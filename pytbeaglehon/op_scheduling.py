#! /usr/bin/env python
"""Classes that accumulate collections of partial calculations to perform
"""
from pytbeaglehon.like_calc_environ import PartialLikeWrapper, NONE_HASH, ProbMatWrapper
class TogglePartialScheduler(object):
    def __init__(self, tree_scorer, model):
        self.tree = tree_scorer._tree
        self.tree_scorer = tree_scorer
        self.like_calc_env = tree_scorer._LCE
        self.model = model
        self._queued_prmat_calcs = []
        self._queued_partial_calcs = []
        self._closed = True
    def end_partial_calculations(self):
        if len(self._queued_prmat_calcs) > 0:
            asrv = self.model.asrv
            if asrv is None:
                nc = 1
                rates = (1.0,)
                asrv_hash = NONE_HASH
            else:
                nc = asrv.num_categories
                rates = asrv.rates
                asrv_hash = asrv.state_hash
            eigen_soln = self.model.eigen_soln
            eigen_state_hash = eigen_soln.state_hash
            for pr_mat_wrapper, edge_len in self._queued_prmat_calcs:
                eff_edge_len_list = [r*edge_len for r in rates]
                ProbMatWrapper.calculate_list(pr_mat_wrapper, 
                                              eigen_soln,
                                              asrv,
                                              eff_edge_len_list,
                                              eigen_hash=eigen_state_hash, 
                                              asrv_hash=asrv_hash)
        if len(self._queued_partial_calcs) > 0:
            fundamental_ops = []
            for el in self._queued_partial_calcs:
                dest, op_type, left_data, left_prob, right_data, right_prob = el
                if op_type == PartialLikeWrapper.NO_LEAF_OP:
                    for n, dest_partial in enumerate(dest):
                        fundamental_ops.append((dest_partial, 
                                                None, # write to this rescaler
                                                None, # read from to this rescaler
                                                left_data[n], 
                                                left_prob[n],
                                                right_data[n],
                                                right_prob[n]))
                elif op_type == PartialLikeWrapper.ONE_LEAF_OP:
                    for n, dest_partial in enumerate(dest):
                        fundamental_ops.append((dest_partial, 
                                                None, # write to this rescaler
                                                None, # read from to this rescaler
                                                left_data[n], 
                                                left_prob[n],
                                                right_data,
                                                right_prob[n]))
                else:
                    for n, dest_partial in enumerate(dest):
                        fundamental_ops.append((dest_partial, 
                                                None, # write to this rescaler
                                                None, # read from to this rescaler
                                                left_data, 
                                                left_prob[n],
                                                right_data,
                                                right_prob[n]))

            PartialLikeWrapper.calc_partials_list(self.like_calc_env, fundamental_ops, to_wait_for=self._LCE_last_queued_dest)
        self._closed = True

    def __del__(self):
        assert(self._closed)

    def queue_prmat(self, node):
        mod = self.model
        if mod not in node._LCE_prob_mat_curr:
            plist = node._LCE_prob_mat_scratch[mod]
            node._LCE_prob_mat_curr[mod] = plist
            self._queued_prmat_calcs.append((plist, node.edge_length))
            return True
        return False

    def add_internal_node_to_partial_calc(self, node):
        self._closed = True
        mod = self.model
        c = node.children
        nc = len(c)
        if nc == 2:
            left_child, right_child = c
            children_dirty = self.queue_prmat(left_child)
            children_dirty = self.queue_prmat(right_child) or children_dirty
            self._LCE_last_queued_dest = []
            if mod not in node._LCE_partial_curr:
                # we need to calculate this partial...
                partial_dest = node._LCE_partial_scratch[mod]
                node._LCE_partial_curr[mod] = partial_dest

                if left_child._LCE_is_internal:
                    assert(mod in left_child._LCE_partial_curr) # child should have been staged.
                    if right_child._LCE_is_internal:
                        assert(mod in right_child._LCE_partial_curr) # child should have been staged.
                        t = (partial_dest,
                             PartialLikeWrapper.NO_LEAF_OP,
                             left_child._LCE_partial_curr[mod],
                             left_child._LCE_prob_mat_curr[mod],
                             right_child._LCE_partial_curr[mod],
                             right_child._LCE_prob_mat_curr[mod])
                    else:
                        t = (partial_dest,
                             PartialLikeWrapper.ONE_LEAF_OP,
                             left_child._LCE_partial_curr[mod],
                             left_child._LCE_prob_mat_curr[mod],
                             right_child._LCE_state_codes,
                             right_child._LCE_prob_mat_curr[mod])
                else:
                    if right_child._LCE_is_internal:
                        assert(mod in right_child._LCE_partial_curr) # child should have been staged.
                        t = (partial_dest,
                             PartialLikeWrapper.ONE_LEAF_OP,
                             right_child._LCE_partial_curr[mod],
                             right_child._LCE_prob_mat_curr[mod],
                             left_child._LCE_state_codes,
                             left_child._LCE_prob_mat_curr[mod])
                    else:
                        t = (partial_dest,
                             PartialLikeWrapper.TWO_LEAF_OP,
                             left_child._LCE_state_codes,
                             left_child._LCE_prob_mat_curr[mod],
                             right_child._LCE_state_codes,
                             right_child._LCE_prob_mat_curr[mod])
                self._queued_partial_calcs.append(t)
                self._LCE_last_queued_dest.extend(partial_dest)
            else:
                assert(not children_dirty) # if we are hear, then we failed to flag the ancestor of a changed prmat as diry
        elif nc > 2:
            raise ValueError("Scoring trees with polytomies is not supported (yet).")
        elif nc == 1:
            raise ValueError("Scoring trees with nodes with outdegree = 1 is not supported (yet)")
        else:
            raise ValueError("add_internal_node_to_partial_calc should only be called for internal nodes")

##############################################################################
##  pytbeaglehon phylogenetic likelihood caluclations using beaglelib.
##
##  Copyright 2010 Mark T. Holder
##  All rights reserved.
##
##  See "LICENSE.txt" for terms and conditions of usage.
##
##############################################################################
