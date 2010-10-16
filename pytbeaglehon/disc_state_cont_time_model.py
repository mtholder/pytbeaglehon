#! /usr/bin/env python
# Copyright (c) 2010 by Mark T. Holder
# (see bottom of file)
"""
Wrapper around a discrete-state, continuous time model of character state change.
"""
_EMPTY_SET = set()

class DiscStateContTimeModel(object):
    def __init__(self, **kwargs):
        self._cmodel = kwargs.get('cmodel')
        self._char_type = kwargs.get('char_type')
        self._model_index = kwargs.get('model_index')
        self._calc_env = kwargs.get('calc_env')
        self._asrv = kwargs.get('asrv')
        self._owns_calc_env = False
        self._changed_params = set()
        self.q_mat = None
        self._q_mat_hash = None
        self._prev_asrv_hash = None
        self._asrv_hash = None
        self._total_state_hash = None
        self._last_asrv_rates_hash = None
        self._prev_state_hash = None
        param_list = kwargs.get('param_list')
        if param_list is not None:
            for p in param_list:
                p.add_listener(self.param_changed)

    def param_changed(self, p):
        '''Adds `p` to this list of changed_parameters.'''
        if self._total_state_hash is not None:
            self._prev_state_hash = self._total_state_hash
        self._total_state_hash = None
        self._changed_params.add(p)

    def get_num_states(self):
        if self._char_type is None:
            return None
        return self._char_type.num_states
    num_states = property(get_num_states)

    def q_mat_is_dirty():
        return (self._changed_params != _EMPTY_SET)

    def asrv_is_dirty():
        self._asrh_hash = self.asrv.state_hash()
        return self._asrh_hash != self._prev_asrv_hash

    def calc_q_mat(self):
        raise NotImplemented()

    def get_q_mat(self):
        if q_mat_is_dirty:
            self._calc_q_mat()
        return self._q_mat

    def set_q_mat(self, v):
        self._changed_params.clear()
        self._total_state_hash = None
        if v is None:
            self._changed_params.add(None) # having a non-empty set assures that the q_mat will be recognized as dirty    
            self._q_mat = None
        else:
            self._q_mat = tuple([tuple([float(i) for i in row]) for row in v])
    q_mat = property(get_q_mat, set_q_mat)
    
    def get_q_mat_hash(self):
        if q_mat_is_dirty():
            self._q_mat_hash = hash(self.q_mat)
        return self._q_mat_hash
    q_mat_hash = property(get_q_mat_hash)

    def state_hash(self):
        if  (self._total_state_hash is None) or  q_mat_is_dirty() or self.asrv_is_dirty():
            self._total_state_hash = hash((id(self), self.q_mat_hash, self.asrv_hash))
            self._prev_asrv_hash = self._asrh_hash
        return self._total_state_hash
        
    def calc_prob_matrices(self, edge_len, eigen_soln_caching=None, prob_mat_caching=None):
        self._incarnate()
        return self._calc_prob_mat(self, edge_len, eigen_soln_caching=eigen_soln_caching, prob_mat_caching=prob_mat_caching)
        

    def _calc_prob_matrices(self, edge_len, eigen_soln_caching=None, prob_mat_caching=None):
        state_id = self.state_hash()
        self._eigen_soln_index = self.calc_env.calc_eigen_soln(q_mat=self.q_mat, state_id=state_id, eigen_soln_caching=eigen_soln_caching)
        self._prob_mat_indices = self._calc_prob_from_eigen(edge_len, self.asrv, state_id=state_id, prob_mat_caching=prob_mat_caching)
        return self._calc_env._get_prob_matrices(self._prob_mat_indices, state_id=state_id)


    def _incarnate(self):
        'Assures that there is a LikeCalcEnvironment associated with this object'
        if (self._cmodel is None) or (self._model_index is None) or (self._calc_env is None):
            if self.num_states is None:
                raise ValueError("DiscStateContTimeModel.num_states must be set before calculations can be performed")
            from pytbeaglehon.like_calc_environ import LikeCalcEnvironment
            self._calc_env = LikeCalcEnvironment(num_leaves=2,
                                                 num_patterns=1,
                                                 num_states=self.num_states,
                                                 num_state_code_arrays=2,
                                                 num_partials=1,
                                                 model_list=[self],
                                                 num_prob_matrices=1,
                                                 num_eigen_storage_structs=1,
                                                 num_rescalings_multipliers=0,
                                                 resource_index=-1)
            self._owns_calc_env = True
    def _reassign_environ(self, calc_env, model_index, cmodel, asrv=None):
        '''Associates the model instance with a new LikeCalcEnvironment

        `model_index` is the new index of this model in that environment
        `cmodel` is a reference to the C object that represents the model.
        '''
        if self._owns_calc_env:
            self._calc_env.release_resources()
            self._owns_calc_env = False
        self._calc_env = calc_env
        self._model_index = model_index
        self._cmodel = cmodel
        self._asrv = asrv
    def get_num_eigen_solutions(self):
        return 1

    def get_num_categories(self):
        a = self.asrv
        if a is None:
            return 1
        return a.num_categories
    num_rate_categories = property(get_num_categories)

    def get_all_submodels(self):
        return [self]
    submodels = property(get_all_submodels)        

    def get_num_prob_models(self):
        return 1
    num_prob_models = property(get_num_prob_models)

class RevDiscStateContTimeModel(DiscStateContTimeModel):

    def __init__(self, **kwargs):
        """Initializes a model. 
        kwargs:
        `r_mat` must be a square matrix (list of lists of floats).  Diagonal
            elements will be ignored. The columns of this matrix will be 
            multiplied by their corresponding equilibirium frequencies.
        `r_upper` can be used instead of `r_mat.` It should correspond to the
            upper triangle of the R-matrix.
        `equil_freq` must be None a list of the equilibrium frequencies.
        """
        r_upper = kwargs.get('r_upper')
        r_mat = kwargs.get('r_mat')
        sf = kwargs.get('state_freq')
        if r_upper:
            if r_mat:
                raise ValueError("r_mat or r_upper cannot both be specified")
            r_mat = _r_upper_to_r_mat(r_upper) 
        elif not r_mat:
            raise ValueError("Either r_mat or r_upper must be given")

        if r_mat:
            self._verify_r_mat(r_mat)
        else:
            raise ValueError("Either r_mat or r_upper must be specified")

        priv_mat = []
        param_set = set()
        for row in r_mat:
            priv_row = list(row)
            for n, cell in enumerate(row):
                try:
                    param_set.update(cell.parameters())
                except:
                    priv_row[n] = MutableFloatParameter(cell)
                
            priv_mat.append(tuple(priv_row))
        self._r_mat = tuple(priv_mat)
        if sf is not None:
            priv_sf = list(sf)
            for cell in sf:
                try:
                    param_set.update(cell.parameters())
                except:
                    priv_sf.append(MutableFloatParameter(cell))
            self.state_freq = tuple(priv_sf)
        else:
            self.state_freq = None

        DiscStateContTimeModel.__init__(self, param_list=param_set, **kwargs)
        if sf is None:
            self.state_freq = tuple([1.0/self.num_states for i in xrange(self.num_states)])

    def get_q_mat(self):
        "Returns a deep copy of the values in the Q-matrix."
        q = []
        for row in self._q_mat:
            q.append(list(row))
        return q

    def get_r_mat(self):
        "Accessor for R-matrix."
        return self._r_mat
    r_mat = property(get_r_mat)
    def get_r_upper(self):
        "Accessor for R-matrix."
        return _r_mat_to_r_upper(self._r_mat)
    r_upper = property(get_r_upper)
    
    def get_state_freq(self):
        "Accessor for state frequencies."
        return self._state_freq
    state_freq = property(get_state_freq)

        
    
    def _recalc_q_mat(self):
        """Uses self._state_freq and self._r_mat to refresh self._q_mat.
        
        As required by the current version of cdsctm_set_q_mat, the Q-matrix
            is rescaled such that each row sums to 0.0 and the weighted 
            average of the diagonal elements is -1.0 (with the weight being the
            associated equilibrium frequency)."""
        _LOG.debug("in self._recalc_q_mat")
        w_mat_sum = 0.0
        for i, rows in enumerate(izip(self._r_mat, self._q_mat)):
            r_row, q_row = rows
            row_sum = 0.0
            for j, to_state in enumerate(izip(r_row, self._state_freq)):
                r_val, stf = [float(p) for p in to_state]
                if i != j:
                    v = r_val * stf
                    assert(r_val >= 0.0)
                    assert(stf >= 0.0)
                    q_row[j] = v
                    row_sum += v
            q_row[i] = -row_sum
            w_mat_sum += row_sum*self._state_freq[i]
        for q_row in self._q_mat:
            for i in xrange(len(q_row)):
                q_row[i] /= w_mat_sum
        _LOG.debug("calling cdsctm_set_q_mat...")
        cdsctm_set_q_mat(self._dsctm, self._q_mat)
        _LOG.debug("back from cdsctm_set_q_mat...")
    
    def _verify_eq_freq_len(self, equil_freq):
        "Raises a ValueError if `equil_freq` is not the correct length"
        nef = len(equil_freq)
        if nef != self._num_states:
            if nef != (1 + self._num_states):
                raise ValueError("The number of states in the r_mat and "\
                             "equil_freq must agree.")

    def _verify_r_mat(self, r_mat):
        """Raises a ValueError if `r_mat` is not the correct shape or has 
        negative values (other than on the diagonal).
        """
        for i, row in enumerate(r_mat):
            if len(row) != self.num_states:
                raise ValueError("The R-matrix must be square.")
            for j, val in enumerate(row):
                if i != j:
                    if val < 0.0:
                        raise ValueError("Off-diagonal elements of the "\
                                         "R-matrix cannot be negative.")
                    if abs(val - float(r_mat[j][i])) > 1.0e-6:
                        raise ValueError("R-matrix must be symmetric")
    def get_num_states(self):
        "Returns the number of states"
        return self._num_states

    def get_adjusted_brlens_list(self, b):
        return [b]

    def get_model_list(self):
        return [self]

    def get_model_probs(self):
        return (1.0,)

    def get_parameters(self):
        return self._parameters

    model_probs = property(get_model_probs)
    r_mat = property(get_r_mat)
    r_upper = property(get_r_upper)
    state_freq = property(get_state_freq)




# Specific models below here...
from pytbeaglehon.disc_char_type import *
class JukesCantorModel(RevDiscStateContTimeModel):
    def __init__(self):
        dna = DNAType()
        _LOG.debug("Created dna type in JukesCantorModel")
        RevDiscreteModel.__init__(self, r_upper=[[1.0, 1.0, 1.0], [1.0, 1.0], [1.0]], char_type=dna, params=[])
        _LOG.debug("called RevDiscreteModel.__init__")


def _r_upper_to_r_mat(r_upper):
    """Convert the upper triangle of a symmetric matrix to the full matrix with 
    0.0 on the diagonal.
    """
    len_upper = len(r_upper)
    r_mat = []
    blank = [0.0]
    for i in xrange(len_upper):
        next_row = blank + r_upper[i]
        blank.append(0.0)
        for j in xrange(i):
            next_row[j] = r_upper[j][i-j-1]
        r_mat.append(next_row)
    for j in xrange(len_upper):
        blank[j] = r_upper[j][len_upper-j-1]
    r_mat.append(blank)
    return r_mat

def _r_mat_to_r_upper(r_mat):
    """Returns the upper triangle (list of lists) of a square matrix `r_mat`"""
    r_upper = []
    for i in xrange(len(r_mat)):
        r_upper.append(r_mat[i][i+1:])
    return r_upper


##############################################################################
##  pytbeaglehon phylogenetic likelihood caluclations using beaglelib.
##
##  Copyright 2010 Mark T. Holder
##  All rights reserved.
##
##  See "LICENSE.txt" for terms and conditions of usage.
##
##############################################################################
