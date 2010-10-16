#! /usr/bin/env python
# Copyright (c) 2005-7 by Mark T. Holder,  University of Kansas
# (see bottom of file)
"""Classes for encapsulating datatypes, and factory for encapsulating data 
structures for characters at the leaves of a tree and the conditional 
likelihood array for a tree.
"""

_DNA_TYPE = None
_AA_TYPE = None
_AA_NO_STOP_TYPE = None

def _create_DNA_type():
    ambig_codes = (("R", "AG"),
                   ("Y","CT"),
                   ("M","AC"),
                   ("W","AT"),
                   ("S","CG"),
                   ("K","GT"),
                   ("V","ACG"),
                   ("H","ACT"),
                   ("D","AGT"),
                   ("B","CGT"),)
    aliases = (("-", "N"),("X", "N"), ("?", "N"))
    return DiscreteCharType("ACGT", 
                                 missing="N", 
                                 ambig_codes=ambig_codes,
                                 aliases=aliases)
def _create_AA_type():
    ambig_codes = (("B", "DN"),
                   ("Z","EQ"), )
    aliases = (("-", "X"), ("?", "X"))
    return DiscreteCharType("ACDEFGHIKLMNPQRSTVWY*", 
                                 missing="X", 
                                 ambig_codes=ambig_codes,
                                 aliases=aliases)
def _create_AA_no_stop_type():
    ambig_codes = (("B", "DN"),
                   ("Z","EQ"), )
    aliases = (("-", "X"), ("?", "X"))
    return DiscreteCharType("ACDEFGHIKLMNPQRSTVWY", 
                                 missing="X", 
                                 ambig_codes=ambig_codes,
                                 aliases=aliases)
def AANoStopType():
    global _AA_NO_STOP_TYPE
    if _AA_NO_STOP_TYPE is None:
        _AA_NO_STOP_TYPE = _create_AA_no_stop_type()
    return _AA_NO_STOP_TYPE

def AAType():
    global _AA_TYPE
    if _AA_TYPE is None:
        _AA_TYPE = _create_AA_type()
    return _AA_TYPE

def DNAType():
    global _DNA_TYPE
    if _DNA_TYPE is None:
        _DNA_TYPE = _create_DNA_type()
    return _DNA_TYPE

class DiscreteCharType(object):
    '''Represents a mapping between a datatype expressed in terms of strings
    and a set of character state indices used internally for likelihood calculations.
    
    Indices [0, nstates) represent the fundamental states. Index nstates is 
    used to indicate complete ambiguity (or missing data).
    
    attributes:
       - `num_states` the number of "fundamental" states
       - `states` a tuple of the labels for the states
       - `all_symbols` is a tuple of all recognize symbols (will not contain both
            cases even if the object ignores case).
      
       - `ignore_case` a bool which is True if case does not matter for state 
            labels.
       - `symbol_to_ind` is a dict mapping a label to the index that the symbol
            represents (note that this does not map ambiguity strings such as 
            {AC} to an index).
       - `state_code_lookup` is maps a state index to a list of the fundamental
            states that are represented by the state.
       - `partial_ambiguity_indices` is a list of indices that map to multiple (but 
            not all of the states)
    methods:
       - `to_indices` is the primary translation function.
       - `has_partial_ambiguity` checks a list of indices for any element of
            partial_ambiguity_indices and returns True if any are found.
    '''

    def __init__(self, states, ambig_codes=(), aliases=None, ignore_case=True, missing=None):
        """Creates an immutable DiscreteCharType object from the list of
        `states`.

        The states list determines the ordinate of states.  `states` must be an
        iterable collection of strings.
        `ambig_codes` must be a iterable collection of key-value pairs. Where the
            keys are new state symbols or state labels, and the values or 
            lists of previously declared states.
        `aliases` is an iterable collection of pairs of strings:
            (alias, real-name).
        If `missing` is specified, it is treated as the first ambig_code
        """
        self._num_states = len(states)
        self._states = tuple(states)
        self._ignore_case = ignore_case
        if self._num_states < 2:
            raise ValueError("The number of states must be greater than 1")
        labels = []
        expansions = []
        symbols_to_state_sets = {}
        symbol_to_ind = {}
        for n, s in enumerate(states):
            if ignore_case:
                s = s.upper()
            if s in symbols_to_state_sets:
                m = "The state %s occurred twice is the state list" % s
                raise ValueError(m)
            
            symbols_to_state_sets[s] = (n,)
            symbol_to_ind[s] = len(labels)
            labels.append(s)
            expansions.append((n,))
        alias_list = aliases and list(aliases) or []
        if missing is not None:
            if missing in symbols_to_state_sets:
                m = "The missing state label (%s) was declared as a state" % missing
                raise ValueError(m)
            ac = [(missing, states)]
            ac.extend(ambig_codes)
        else:
            ac = ambig_codes
        self._partial_ambig_indices = []
        for key, val in ac:
            if ignore_case:
                key = key.upper()
            if key in symbols_to_state_sets:
                m = "The ambiguity code (%s) was declared as a state" % key
                raise ValueError(m)
            expanded = set()
            for v_el in val:
                if ignore_case:
                    v_el = v_el.upper()
                try:
                    e_el = symbols_to_state_sets[v_el]
                    for e in e_el:
                        expanded.add(e)
                except KeyError:
                    raise ValueError("The ambiguity code expansion (for code %s) contains an unknown label (%s)" % (key, v_el))
            exp_list = list(expanded)
            exp_list.sort()
            texp = tuple(exp_list)
            if texp in expansions:
                i = expansions.index(texp)
                alias_list.append((key, labels[i]))
            else:
                symbols_to_state_sets[key] = texp
                symbol_to_ind[key] = len(labels)
                if len(texp) > 1 and len(texp) < self.num_states:
                    self._partial_ambig_indices.append(len(labels))
                labels.append(key)
                expansions.append(texp)
        for key, val in alias_list:
            if ignore_case:
                key = key.upper()
                val = val.upper()
            if key in symbols_to_state_sets:
                m = "The alias (%s) was declared as a state" % key
                raise ValueError(m)
            if val not in symbols_to_state_sets:
                m = "The alias %s maps to an unknown state (%s)" % (key, val)
                raise ValueError(m)
            symbol_to_ind[key] = symbol_to_ind[val]
        self._symbol_to_ind = symbol_to_ind
        self._all_symbols = tuple(labels)
        self._state_code_lookup = tuple(expansions)
        expansions = []
        for i in self._state_code_lookup:
            expansions.append(list(i))
        self._partial_ambig_indices = tuple(self._partial_ambig_indices)

    def to_indices(self, seq):
        "Converts a sequence of character symbols to the corresponding indices."
        l = []
        sti = self._symbol_to_ind
        for s in seq:
            l.append(sti[s])
        return l

    def get_num_states(self):
        return self._num_states
    num_states = property(get_num_states)

    def get_states(self):
        return self._states
    states = property(get_states)

    def get_ignore_case(self):
        return self._ignore_case
    ignore_case = property(get_ignore_case)

    def get_symbol_to_ind(self):
        return dict(self._symbol_to_ind)
    symbol_to_ind = property(get_symbol_to_ind)

    def get_all_symbols(self):
        return self._all_symbols
    all_symbols = property(get_all_symbols)

    def get_state_code_lookup(self):
        return self._state_code_lookup
    state_code_lookup = property(get_state_code_lookup)

    def get_partial_ambiguity_indices(self):
        return self._partial_ambig_indices
    partial_ambiguity_indices = property(get_partial_ambiguity_indices)

    def has_partial_ambiguity(self, indices):
        pai = set(self.partial_ambiguity_indices)
        for el in indices:
            if el in pai:
                return True
        return False
        
    

##############################################################################
##  pytbeaglehon phylogenetic likelihood caluclations using beaglelib.
##
##  Copyright 2010 Mark T. Holder
##  All rights reserved.
##
##  See "LICENSE.txt" for terms and conditions of usage.
##
##############################################################################
