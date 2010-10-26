#! /usr/bin/env python
"""Wrapper class for scoring a tree. Generally obtained from a LikeCalcEnvironment 
instance through the tree_scorer method.

Expectations of a tree object:
    
    postorder_internal_node_iter returns node objects for all internal nodes 
            in postorder traversal order.
    
    TreeScorer initialization will call postorder_internal_node_iter and this 
        should sweep over *all* internal nodes that will be encountered (and the
        children of these nodes will include all leaves).
        
    Nodes do not change from internal to leaf or vice versa.
    
    
Expectations of a node objects:
    "children" attribute (or property) is a list of two children Nodes or an
            empty list
    
    "parent" attribute (or property) should be another Node object (or None 
            for the root)
    
    No attributes used by client code start with _LCE_
    
    leaves have integer "leaf_index" which is in [0, num_leaves) and agrees with
        the order that data was entered in the LCE.

    every node should have a edge_length attribute (or property) with a non-negative
        float.
"""

class TreeScorer(object):
    def __init__(self, like_calc_env, tree):
        self._tree = tree
        self._LCE = like_calc_env
        self._edge_to_edge_wrapper = {}
        self._edge_to_wrapper = {}
        self._node_to_wrapper = {}
        self._changed_edges = set()
        self._changed_nodes = set()
        self._entire_tree_dirty = True
        self._changed_models = set(like_calc_env.model_list)
        self._num_models = len(like_calc_env.model_list)
        self._cached_model_to_score = {}
        self.initialize_tree()
    def initialize_tree(self):
        pass
        
        
        
    def get_entire_tree_is_dirty(self):
        return self._entire_tree_dirty
    entire_tree_is_dirty = property(get_entire_tree_is_dirty)
    def __call__(self):
        if self.entire_tree_is_dirty:
            self._changed_models = set(self._LCE.model_list)
        changed_model_list = list(self._changed_models)
        num_mods_to_update = len(changed_model_list)
        
        for cm in changed_model_list:
            self._cached_model_to_score[cm] = self._calc_full_traversal_lnL_for_model(cm)
            self._changed_models.remove(cm)

        if num_mods_to_update < self._num_models:
            raise NotImplementedError("updates of part of the tree are not supported yet.")

        self._entire_tree_dirty = False
        self._changed_edges.clear()
        self._changed_nodes.clear()

        self._lnL = sum([i for i in self._cached_model_to_score.itervalues()])
        return self._lnL

    def _calc_full_traversal_lnL_for_model(self, model):
        self._LCE.start_partial_calculations(model)
        try:
            for nd in self._tree.postorder_internal_node_iter():
                 self._LCE.add_internal_node_to_partial_calc(model, nd)
        finally:
            self._LCE.end_partial_calculations(model)
        return 0.0

class TogglePartialTreeScorer(object):
    '''Assumes that there are:
            - enough partials for every internal node to have two
            - enough prob mats for every edge to have two.
       This enables an efficient "propose" followed by "accept" or "reject"
       API.

       Note that when reject is called, the client must then call:
            1. start_revert()
            2. parameter changing moves to return the tree to the previous state,
            3. end_revert()
       So that the changes to restore the tree do not result in "dirty-ing" of 
        partials.
    '''

    def initialize_tree(self):
        LCE = self._LCE
        tips = []
        for nd in self._tree.postorder_internal_node_iter():
            if nd.parent is None:
                nd._LCE_prob_mat_stored = None
                nd._LCE_prob_mat_scratch = None
            else:
                nd._LCE_prob_mat_stored = LCE._prob_mat_cache.get_writable_object()
                nd._LCE_prob_mat_scratch = LCE._prob_mat_cache.get_writable_object()
                nd._LCE_edge_len_stored = nd.edge_length
                nd._LCE_edge_len_scratch = nd.edge_length
            nd._LCE_is_internal = True
            nd._LCE_partial_stored = LCE._partial_cache.get_writable_object()
            nd._LCE_partial_scratch = LCE._partial_cache.get_writable_object()
            for c in nd.children:
                if not bool(c.children):
                    tips.append(c)
        leaf_inds = set()
        for nd in tips:
            if nd.leaf_index in leaf_inds:
                raise ValueError("The leaf_index %s repeated!" % str(nd.leaf_index))
            leaf_inds.add(nd.leaf_index)
            nd._LCE_prob_mat_stored = LCE._prob_mat_cache.get_writable_object()
            nd._LCE_prob_mat_scratch = LCE._prob_mat_cache.get_writable_object()
            nd._LCE_edge_len_stored = nd.edge_length
            nd._LCE_edge_len_scratch = nd.edge_length
            nd._LCE_is_internal = False
            nd._LCE_buffer_index = nd.leaf_index
            nd._LCE_state_codes = LCE._wrap_state_code_array[nd.leaf_index]
            assert(nd._LCE_state_codes.is_calculated)
        tips.sort(cmp=lambda x, y: cmp(x.leaf_index, y.leaf_index))
        max_leaf_ind = tips[-1].leaf_index
        self._tips = [None] *(1+max_leaf_ind)
        for nd in tips:
            self._tips[nd.leaf_index] = nd
            
    
##############################################################################
##  pytbeaglehon phylogenetic likelihood caluclations using beaglelib.
##
##  Copyright 2010 Mark T. Holder
##  All rights reserved.
##
##  See "LICENSE.txt" for terms and conditions of usage.
##
##############################################################################
