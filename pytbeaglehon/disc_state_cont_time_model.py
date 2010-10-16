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
        self._num_states = kwargs.get('num_states')
        self._char_type = kwargs.get('char_type')
        if (self._char_type is not None):
            if (self._num_states is not None) and (self._num_states != self._char_type.num_states):
                raise ValueError("If 'num_states' and 'char_type' are both specified, then the number of states must agree")
            self._num_states = self._char_type.num_states
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

class RevDiscStateContTimeModel(DiscStateContTimeModel):
    pass

class JukesCantorModel(RevDiscStateContTimeModel):
    def __init__(self):
        dna = DNAType()
        _LOG.debug("Created dna type in JukesCantorModel")
        RevDiscreteModel.__init__(self, r_upper=[[1.0, 1.0, 1.0], [1.0, 1.0], [1.0]], char_type=dna, params=[])
        _LOG.debug("called RevDiscreteModel.__init__")

##############################################################################
##  pytbeaglehon phylogenetic likelihood caluclations using beaglelib.
##
##  Copyright 2010 Mark T. Holder
##  All rights reserved.
##
##  See "LICENSE.txt" for terms and conditions of usage.
##
##############################################################################
