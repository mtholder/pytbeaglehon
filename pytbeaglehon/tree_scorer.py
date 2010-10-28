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
        self.model_list = like_calc_env.model_list
        self.initialize_tree()

    def initialize_tree(self):
        pass
        
        
        
    def get_entire_tree_is_dirty(self):
        return self._entire_tree_dirty
    entire_tree_is_dirty = property(get_entire_tree_is_dirty)

    def __call__(self):
        return get_ln_L(self)

    def get_ln_L(self):
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
class TogglePartialTreeScorer(TreeScorer):
    '''Assumes that there are:
            - enough partials for every internal node to have two
            - enough prob mats for every edge to have two.
       This enables an efficient API of repeatedly:
            1. propose-new-state (by calling ...param_changed methods) followed by,
            2. "accept" or "reject" call

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
                nd._LCE_edge_len_stored = None
                nd._LCE_edge_len_scratch = None
            else:
                nd._LCE_prob_mat_stored = {}
                nd._LCE_prob_mat_scratch = {}
                for mod in self.model_list:
                    nd._LCE_prob_mat_stored[mod] = [LCE._prob_mat_cache.get_writable_object() for i in range(mod.num_rate_categories)]
                    nd._LCE_prob_mat_scratch[mod] = [LCE._prob_mat_cache.get_writable_object() for i in range(mod.num_rate_categories)]
                nd._LCE_edge_len_stored = nd.edge_length
                nd._LCE_edge_len_scratch = nd.edge_length
            nd._LCE_edge_len_curr = None
            nd._LCE_prob_mat_curr = {}
            nd._LCE_partial_curr = {}
            nd._LCE_is_internal = True
            nd._LCE_partial_stored = {}
            nd._LCE_partial_scratch = {}
            for mod in self.model_list:
                nd._LCE_partial_stored[mod] = [LCE._partial_cache.get_writable_object() for i in range(mod.num_rate_categories)]
                nd._LCE_partial_scratch[mod] = [LCE._partial_cache.get_writable_object() for i in range(mod.num_rate_categories)]
            for c in nd.children:
                if not bool(c.children):
                    tips.append(c)
        leaf_inds = set()
        for nd in tips:
            if nd.leaf_index in leaf_inds:
                raise ValueError("The leaf_index %s repeated!" % str(nd.leaf_index))
            leaf_inds.add(nd.leaf_index)
            nd._LCE_prob_mat_stored = {}
            nd._LCE_prob_mat_scratch = {}
            for mod in self.model_list:
                nd._LCE_prob_mat_stored[mod] = [LCE._prob_mat_cache.get_writable_object() for i in range(mod.num_rate_categories)]
                nd._LCE_prob_mat_scratch[mod] = [LCE._prob_mat_cache.get_writable_object() for i in range(mod.num_rate_categories)]
            nd._LCE_prob_mat_curr = {}
            nd._LCE_edge_len_stored = nd.edge_length
            nd._LCE_edge_len_scratch = nd.edge_length
            nd._LCE_edge_len_curr = None
            nd._LCE_is_internal = False
            nd._LCE_buffer_index = nd.leaf_index
            nd._LCE_state_codes = LCE._wrap_state_code_array[nd.leaf_index]
            if not nd._LCE_state_codes.is_calculated:
                raise ValueError("State codes (data for the leaf of a tree) has not been specified for leaf_index=%d" % nd.leaf_index)
        tips.sort(cmp=lambda x, y: cmp(x.leaf_index, y.leaf_index))
        max_leaf_ind = tips[-1].leaf_index
        self._tips = [None] *(1+max_leaf_ind)
        for nd in tips:
            self._tips[nd.leaf_index] = nd
        self._scheduler = None
        self.get_ln_L()
        self.accept()

    
    def accept(self):
        for nd in self._tree.postorder_internal_node_iter():
            if nd._LCE_edge_len_curr is nd._LCE_edge_len_scratch:
                nd._LCE_edge_len_stored, nd._LCE_edge_len_scratch = nd._LCE_edge_len_scratch, nd._LCE_edge_len_stored
            for mod in self.model_list:
                if nd._LCE_prob_mat_curr[mod] is nd._LCE_prob_mat_scratch[mod]:
                    nd._LCE_prob_mat_stored[mod], nd._LCE_prob_mat_scratch[mod] = nd._LCE_prob_mat_scratch[mod], nd._LCE_prob_mat_stored[mod]
                if nd._LCE_partial_curr[mod] is nd._LCE_partial_scratch[mod]:
                    nd._LCE_partial_stored[mod], nd._LCE_partial_scratch[mod] = nd._LCE_partial_scratch[mod], nd._LCE_partial_stored[mod]

    def revert(self):
        for nd in self._tree.postorder_internal_node_iter():
            nd._LCE_edge_len_curr = nd._LCE_edge_len_stored
            nd._LCE_prob_mat_curr = dict(nd._LCE_prob_mat_stored)
            nd._LCE_partial_curr = dict(nd._LCE_partial_stored)

    def _calc_full_traversal_lnL_for_model(self, model):
        assert(self._scheduler is None)
        from pytbeaglehon.op_scheduling import TogglePartialScheduler
        self._scheduler = TogglePartialScheduler(self, model)
        try:
            for node in self._tree.postorder_internal_node_iter():
                 self._scheduler.add_internal_node_to_partial_calc(node)
        finally:
            self._scheduler.end_partial_calculations()
            self._scheduler = None

        return 0.0


def create_toggle_partial_tree_scorer(model_list, 
                                      data,
                                      tree,
                                      num_extra_nodes=0
                                      ):
    '''Simple constructor for a TogglePartialTreeScorer and LikeCalcEnvironment
    to support it.
    
    Infers:
        data_type (num_states) and asrv from the models in model_list
        num_leaves, num_state_code_arrays, num_partials for len(data)
        num_patterns from len(data[0]) 
        
        
        
        Assumes that you will want a enough prob_matrix for every edge in a rooted,
            binary tree to have two sets of matrices for each model/rate-category 
            combination.
        Assumes that you want only two eigen solution per model.

        Assumes that you want only two partials for internal nodes solution per
            model and `num_extra_partials` in addition to this (thus if you want
            one "extra" node that can be swapped in and out of the tree you 
            should use num_extra_partials=1, and then manually load up the 
            _LCE_xxx attributes for that node.

        Assumes that you want one rescaling array for every 6 edges (every 4 leaves)
    '''
    from pytbeaglehon.like_calc_environ import LikeCalcEnvironment
    asrv_list = []
    num_model_rate_cats = 0
    num_leaves = len(data)
    num_patterns = len(data[0])
    num_models = len(model_list) #TODO more generic for mixtures!
    for model in model_list:
        a = model.asrv
        if a is None:
            num_model_rate_cats += 1
        else:
            num_model_rate_cats += a.num_categories
    num_internals = (num_leaves - 1)
    num_nodes = num_internals + num_leaves
    LCE = LikeCalcEnvironment(model_list=model_list,
                               num_patterns=num_patterns,
                               num_leaves=num_leaves,
                               num_state_code_arrays=num_leaves,
                               num_partials=(num_internals + num_extra_nodes)*2*num_model_rate_cats,
                               num_prob_matrices=(num_nodes + num_extra_nodes)*2*num_model_rate_cats,
                               num_eigen_storage_structs=2*num_models,
                               num_rescalings_multipliers= 2*(1 + num_leaves//4))
    for n, row in enumerate(data):
        LCE.set_state_code_array(n, row)
    scorer = TogglePartialTreeScorer(LCE, tree)
    return scorer

##############################################################################
##  pytbeaglehon phylogenetic likelihood caluclations using beaglelib.
##
##  Copyright 2010 Mark T. Holder
##  All rights reserved.
##
##  See "LICENSE.txt" for terms and conditions of usage.
##
##############################################################################
