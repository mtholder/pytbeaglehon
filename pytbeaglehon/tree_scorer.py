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
    "children" attribute returns a list of two children or an empty list
    
    No attributes used by client code start with _LCE_
    
    leaves have integer "leaf_index" which is in [0, num_leaves) and agrees with
        the order that data was entered in the LCE.
    
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
        
        tips = []
        for nd in self._tree.postorder_internal_node_iter():
            nd._LCE_prmat_index = None
            nd._LCE_edge_state_cache = None
            nd._LCE_is_internal = True
            nd._LCE_save_partials = True
            nd._LCE_buffer_index = None
            nd._LCE_buffer_state_id = None
            nd._LCE_partial_state_cache = None
            nd._LCE_edge_state_cache = None
            for c in nd.children:
                if not bool(c.children):
                    tips.append(c)
        leaf_inds = set()
        for nd in tips:
            if nd.leaf_index in leaf_inds:
                raise ValueError("The leaf_index %s repeated!" % str(nd.leaf_index))
            leaf_inds.add(nd.leaf_index)
            nd._LCE_prmat_index = None
            nd._LCE_edge_state_cache = None
            nd._LCE_is_internal = False
            nd._LCE_buffer_index = nd.leaf_index
            nd._LCE_buffer_state_id = str('leaf%d' % nd.leaf_index)
            nd._LCE_edge_state_cache = None
        tips.sort(cmp=lambda x, y: cmp(x.leaf_index, y.leaf_index))
        max_leaf_ind = tips[-1].leaf_index
        self._tips = [None] *(1+max_leaf_ind)
        for nd in tips:
            self._tips[nd.leaf_index] = nd
            
        
        
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
##############################################################################
##  pytbeaglehon phylogenetic likelihood caluclations using beaglelib.
##
##  Copyright 2010 Mark T. Holder
##  All rights reserved.
##
##  See "LICENSE.txt" for terms and conditions of usage.
##
##############################################################################
