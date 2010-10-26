#! /usr/bin/env python
from itertools import izip
from pytbeaglehon.ccore.disc_state_cont_time_model import \
    cpytbeaglehon_init, cpytbeaglehon_free, cget_num_comp_resources, \
    cget_comp_resource_info, cget_model_list, cdsctm_set_q_mat, \
    cdsctm_calc_eigens, cdsctm_calc_pr_mats, cdsctm_get_pr_mats, \
    cdsctm_set_state_code
from pytbeaglehon import DiscStateContTimeModel
from pytbeaglehon import get_logger, CachingFacets
_LOG = get_logger(__name__)
_EMPTY_SET = set()
_EMPTY_DICT = {}


class BeagleResourceFlags:
    PRECISION_SINGLE    = 1 << 0 #  Single precision computation 
    PRECISION_DOUBLE    = 1 << 1 #  Double precision computation 

    COMPUTATION_SYNCH   = 1 << 2 #  Synchronous computation (blocking) 
    COMPUTATION_ASYNCH  = 1 << 3 #  Asynchronous computation (non-blocking) 
    
    EIGEN_REAL          = 1 << 4 #  Real eigenvalue computation 
    EIGEN_COMPLEX       = 1 << 5 #  Complex eigenvalue computation 

    SCALING_MANUAL      = 1 << 6 #  Manual scaling 
    SCALING_AUTO        = 1 << 7 #  Auto-scaling on 
    SCALING_ALWAYS      = 1 << 8 #  Scale at every updatePartials 
    SCALING_DYNAMIC     = 1 << 19 #  Manual scaling with dynamic checking  
    
    SCALERS_RAW         = 1 << 9 #  Save raw scalers 
    SCALERS_LOG         = 1 << 10 #  Save log scalers 
    
    VECTOR_SSE          = 1 << 11 #  SSE computation 
    VECTOR_NONE         = 1 << 12 #  No vector computation 
    
    THREADING_OPENMP    = 1 << 13 #  OpenMP threading 
    THREADING_NONE      = 1 << 14 #  No threading 
    
    PROCESSOR_CPU       = 1 << 15 #  Use CPU as main processor 
    PROCESSOR_GPU       = 1 << 16 #  Use GPU as main processor 
    PROCESSOR_FPGA      = 1 << 17 #  Use FPGA as main processor 
    PROCESSOR_CELL      = 1 << 18 #  Use Cell as main processor 

# add to_flag_name and to_flag_number dictionaries to BeagleResourceFlags
_name_to_num = {}
_num_to_name = {}
for _k, _v in BeagleResourceFlags.__dict__.items():
    if isinstance(_k, str) and (isinstance(_v, int) or isinstance(_v, long)):
        _name_to_num[_k] = _v
        _num_to_name[_v] = _k
BeagleResourceFlags.to_flag_name = _num_to_name
BeagleResourceFlags.to_flag_number = _name_to_num
del _num_to_name
del _name_to_num



def minimal_LCE(model_list, data):
    '''Simple constructor for a LikeCalcEnvironment that infers:
        data_type (num_states) and asrv from the models in model_list
        num_leaves, num_state_code_arrays, num_partials for len(data)
        num_patterns from len(data[0])
        
        Assumes that you will want a enough prob_matrix for every edge in a rooted,
            binary tree to have a set of matrices for each model/rate-category 
            combination.
        Assumes that you want only one eigen solution per model.

        Assumes that you want one rescaling array for every 6 edges (every 4 leaves)
    '''
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
        
    LCE = LikeCalcEnvironment(model_list=model_list,
                               num_patterns=num_patterns,
                               num_leaves=num_leaves,
                               num_state_code_arrays=num_leaves,
                               num_partials=(num_leaves - 1)*num_model_rate_cats,
                               num_prob_matrices=num_model_rate_cats,
                               num_eigen_storage_structs=num_models,
                               num_rescalings_multipliers= 1 + num_leaves//4)
    for n, row in enumerate(data):
        LCE.set_state_code_array(n, row)
    return LCE

_NONE_HASH = '<None>'

class BufferWrapper(object):
    '''Base class for Python objects that wrap a beagle buffer.  Holds a
    reference to the LikeCalcEnvironment and the index of the struct within
    that context.
    '''
    def __init__(self, index, like_calc_env):
        self.index = index
        self.like_calc_env = like_calc_env
        self._is_calculated = False
    def set_calculated(self):
        self._is_calculated = True
    def get_is_calculated(self):
        return self._is_calculated
    def clear(self):
        self._is_calculated = False
    is_calculated = property(get_is_calculated)

class EigenSolutionWrapper(BufferWrapper):

    def __init__(self, index, like_calc_env):
        BufferWrapper.__init__(self, index=index, like_calc_env=like_calc_env)
        self._instance_hash_format = 'ES-%d-%d(%%s)' % (id(self), index)
        self.clear()
    def clear(self):
        self._is_calculated = False
        self._model = None 
        self._model_hash = None
        self._state_hash = None
        self._instance_hash = None

    def calc_hash(model_state_hash):
        return '%s' % model_state_hash
    calc_hash = staticmethod(calc_hash)

    def get_state_hash(self):
        if self._state_hash is None:
            if self._model_hash is None:
                raise ValueError('EigenSolutionWrapper with empty model is not hashable')
            self._state_hash = EigenSolutionWrapper.calc_hash(self._model_hash)
        return self._state_hash
    state_hash = property(get_state_hash)

    def get_instance_hash(self):
        if self._instance_hash is None:
            if self._model_hash is None:
                raise ValueError('EigenSolutionWrapper with empty model is not hashable')
            self._instance_hash = self._instance_hash_format % self._model_hash
        return self._state_hash
    instance_hash = property(get_instance_hash)
  
    def calculate(self, model, model_state_hash=None):
        cdsctm_set_q_mat(model.cmodel, model.q_mat)
        _LOG.debug("Calling cdsctm_calc_eigens(%d, %d)" % (id(model), self.index))
        cdsctm_calc_eigens(model.cmodel, self.index)
        self._state_hash = None
        self._model = model
        if model_state_hash is None:
            self._model_hash = model.state_hash
        else:
            assert(model.state_hash == model_state_hash) # doublechecking
            self._model_hash = model_state_hash

class ProbMatWrapper(BufferWrapper):
    _hash_format = '%s-%s-%d-%s'
    def calc_hash(eigen_state_hash, asrv_hash, asrv_categ, edge_len_hash):
        return ProbMatWrapper._hash_format % (eigen_state_hash, 
                                              asrv_hash,
                                              asrv_categ,
                                              edge_len_hash)
    calc_hash = staticmethod(calc_hash)
    def __init__(self, index, like_calc_env):
        BufferWrapper.__init__(self, index=index, like_calc_env=like_calc_env)
        self._instance_hash_format = ('PM-%d-%d' % (id(self), index)) + ProbMatWrapper._hash_format
        self.clear()
    def clear(self):
        self._is_calculated = False
        self._eigen_solution = None
        self._eigen_solution_hash = None
        self._asrv = None
        self._asrv_hash = None
        self._asrv_categ = None
        self._edge_length = None
        self._edge_length_hash = None
        self._state_hash = None
        self._instance_hash = None
    
    def get_state_hash(self):
        if self._state_hash is None:
            if self._eigen_solution_hash is None:
                raise ValueError('ProbMatWrapper without an eigen solution is not hashable')
            if self._asrv_hash is None:
                raise ValueError('ProbMatWrapper without an asrv object is not hashable')
            if self._asrv_categ is None:
                raise ValueError('ProbMatWrapper without an asrv category is not hashable')
            if self._edge_length_hash is None:
                raise ValueError('ProbMatWrapper without an edge_length is not hashable')
            self._state_hash = ProbMatWrapper._hash_format % (self._eigen_solution_hash,
                                                   self._asrv_hash,
                                                   self._asrv_categ,
                                                   self._edge_length_hash)
        return self._state_hash
    state_hash = property(get_state_hash)

    def get_instance_hash(self):
        if self._instance_hash is None:
            if self._eigen_solution_hash is None:
                raise ValueError('ProbMatWrapper without an eigen solution is not hashable')
            if self._asrv_hash is None:
                raise ValueError('ProbMatWrapper without an asrv object is not hashable')
            if self._asrv_categ is None:
                raise ValueError('ProbMatWrapper without an asrv category is not hashable')
            if self._edge_length_hash is None:
                raise ValueError('ProbMatWrapper without an edge_length is not hashable')
            self._instance_hash = self._instance_hash_format % (self._eigen_solution_hash,
                                                   self._asrv_hash,
                                                   self._asrv_categ,
                                                   self._edge_length_hash)
        return self._instance_hash
    instance_hash = property(get_instance_hash)
    def calculate_list(pr_wrap_list, eigen_soln, asrv, eff_edge_len_list, eigen_hash=None, asrv_hash=None, eff_edge_len_hash=None):
        lce = pr_wrap_list[0].like_calc_env
        ind_list = []
        for p in pr_wrap_list:
            assert(p.like_calc_env is lce)
            ind_list.append(p.index)

        _LOG.debug("Calling cdsctm_calc_pr_mats(%d, %d, %s, %s)" % (lce._handle, eigen_soln.index, str(eff_edge_len_list), str(ind_list)))
        cdsctm_calc_pr_mats(lce._handle, eigen_soln.index, eff_edge_len_list, ind_list)
        
        for n, p in enumerate(pr_wrap_list):
            p._eigen_solution = eigen_soln
            p._eigen_solution_hash = eigen_hash
            p._asrv = asrv
            p._asrv_hash = asrv_hash
            p._asrv_categ = n
            p._edge_length = eff_edge_len_list[n]
            p._edge_length_hash = eff_edge_len_hash[n]
            p._state_hash = None
            p._instance_hash = None
        return pr_wrap_list
    calculate_list = staticmethod(calculate_list)


class PartialLikeWrapper(BufferWrapper):
    def __init__(self, index, like_calc_env):
        BufferWrapper.__init__(self, index=index, like_calc_env=like_calc_env)
        self.revision_index = None # stores the number of times that the wrapped object has changed -- but reversions are allowed
        self.next_revision_index = 0 # next unique identifier (will 1+self.revision_index if the current state is not a reversion to a previous state
        self._state_hash_format = 'PL-%d-%d-%%d' % (id(self), index)
        self.full_hash_format = 'PL-%d-%d-%%d(%%s-%%s+%%s-%%s)' % (id(self), index)
        self.clear()

    def clear(self):
        self._is_calculated = False
        self._left_data_hash = None
        self._left_prmat = None
        self._left_prmat_hash = None
        self._right_data_hash = None
        self._right_prmat = None
        self._right_prmat_hash = None
        self._state_hash = None
        self._full_state_hash = None
        self.revision_index = None

    def set_calculated(self):
        self.revision_index = self.next_revision_index
        self.next_revision_index += 1
        
    def get_full_state_hash(self):
        if self._full_state_hash is None:
            if self.revision_index is None:
                raise ValueError('ProbMatWrapper that has not been calculated is not hashable')
            if self._left_data_hash is None:
                raise ValueError('ProbMatWrapper without an left child data is not hashable')
            if self._left_prmat_hash is None:
                raise ValueError('ProbMatWrapper without an left child probability matrix is not hashable')
            if self._right_data_hash is None:
                raise ValueError('ProbMatWrapper without an right child data is not hashable')
            if self._right_prmat_hash is None:
                raise ValueError('ProbMatWrapper without an right child probability matrix is not hashable')
            self._full_state_hash = self.full_hash_format % (self.revision_index,
                                                   self._left_data_hash,
                                                   self._left_prmat_hash,
                                                   self._right_data_hash,
                                                   self._right_prmat_hash)
        return self._full_state_hash
    full_state_hash = property(get_full_state_hash)
    def get_state_hash(self):
        if self._state_hash is None:
            if self.revision_index is None:
                raise ValueError('ProbMatWrapper that has not been calculated is not hashable')
            self._state_hash = self._state_hash_format % (self.revision_index)
        return self._state_hash
    state_hash = property(get_state_hash)
 

class StateCodeArrayWrapper(BufferWrapper):
    def __init__(self, index, like_calc_env):
        BufferWrapper.__init__(self, index=index, like_calc_env=like_calc_env)
        self._leaf_index = index
        self._state_hash_format = 'SC-%d-%d-%%d-%%d' % (id(self), index)
        self.next_revision_index = 0 # next unique identifier (will 1+self.revision_index if the current state is not a reversion to a previous state
        self.clear()

    def clear(self):
        self._is_calculated = False
        self._state_hash = None
        self.revision_index = None # stores the number of times that the wrapped object has changed -- but reversions are allowed

    def set_calculated(self):
        self.revision_index = self.next_revision_index
        self.next_revision_index += 1

    def get_state_hash(self):
        if self._state_hash is None:
            if self.revision_index is None:
                raise ValueError('StateCodeArrayWrapper that has not been set is not hashable')
            self._state_hash = self._state_hash_format % (self._leaf_index, self.revision_index)
        return self._state_hash
    state_hash = property(get_state_hash)
    full_state_hash = property(get_state_hash)


class RescalingMultiplier(BufferWrapper):
    def __init__(self, index, like_calc_env):
        BufferWrapper.__init__(self, index=index, like_calc_env=like_calc_env)
        self.hash_format = 'RM-%d-%d(%%s)' % (id(self), index)
        self.clear()

    def clear(self):
        self._is_calculated = False
        self._state_hash = None
        self._partial_wrapper = None
        self._partial_wrapper_hash = None

    def get_state_hash(self):
        if self._state_hash is None:
            if self._partial_wrapper_hash is None:
                raise ValueError('RescalingMultiplier that has not been calculated is not hashable')
            self._state_hash = self.hash_format % (self._partial_wrapper_hash)
        return self._state_hash
    state_hash = property(get_state_hash)

class CalculatedCache(object):
    
    def __init__(self, wrappers, obj_name):
        self.obj_name = obj_name
        self._free = set(wrappers)
        self._saved = {}
        self._calculated = set()
        self._state_to_wrapper = {}
        self._queued = set()

    def get_from_cache(self, state_hash):
        return self._state_to_wrapper.get(state_hash)

    def save_obj(self, o):
        '''increments the reference count on o'''
        self._queued.discard(o)
        o.set_calculated()
        n = self._saved.setdefault(o, 0)
        self._saved[o] = (n + 1)

    def get_writable_object(self, o=None):
        '''Returns a free object, and "tells" the wrapper to clear itself.
        
        Following this call, the caller must call exactly one of the following:
            - `flag_as_calculated` (to signal the calculation was successful, but
                the object does not need to be stored long term).
            - `save_obj` to save the object, or
            - `release` to return the object (as uncalculated) to the free pool
        '''
        if o is None:
            try:
                o = self._free.pop()
            except KeyError:
                try:
                    o = self._calculated.pop()
                except KeyError:
                    raise ValueError("All %s instances are locked" % self.obj_name)
            o.clear()
            self._queued.add(o)
        else:
            self.make_writable(o)
        return o

    def flag_as_calculated(self, o):
        assert(o in self._queued)
        assert(o not in self._saved)
        o.set_calculated()
        self._queued.discard(o)
        self._calculated.add(o)

    def save_obj(self, o):
        assert(o in self._queued)
        assert(o not in self._saved)
        o.set_calculated()
        self._queued.discard(o)
        self._saved[o] = 1

    def release(self, o):
        assert(o in self._queued)
        assert(o not in self._free)
        self._queued.discard(o)
        self._free.add(o)

    def make_writable(self, o):
        r = self._saved.get(o)
        if r is None:
            if o in self._calculated:
                self._calculated.discard(o)
            else:
                assert(o in self._free)
        elif r == 1:
            del self._saved[o]
        else:
            return self.get_writable_object()
        o.clear()
        self._queued.add(o)
        return o
        
        

class LikeCalcEnvironment(object):
    _CALC_ACTION = 0
    _FREE_ACTION = 1
    
    def get_num_comp_resources():
        return cget_num_comp_resources()
    get_num_comp_resources = staticmethod(get_num_comp_resources)

    def query_comp_resource_info(resourceIndex):
        return cget_comp_resource_info(resourceIndex)
    query_comp_resource_info = staticmethod(query_comp_resource_info)

    def get_comp_resource_info(self):
        if self._resource_index is None:
            raise ValueError("resource_index (the index of the computational resource) must be set before comp_resource_info can be accessed")
        _LOG.debug("calling query_comp_resource_info for %d" % self._resource_index)
        return LikeCalcEnvironment.query_comp_resource_info(self._resource_index)
    comp_resource_info = property(get_comp_resource_info)

    def __init__(self, **kwargs):
        """Creates an new instance of a likelihood calculation context 
        this is necessary because beagle does not support adding memory to an
        instance after initialization.
        
        keyword arguments (these can also be set as attributes of the object
          as long as they are set before the initialization of the beagle
          environment is triggered by requesting the model objects or some
          other calculation that requires the instance be running):
            `num_leaves`
            `num_patterns`
            `pattern_weight_list`
            `num_states`
            `num_state_code_arrays`
            `num_partials`
            `num_inst_rate_matrices`
            `asrv_list`
            `num_prob_matrices`
            `num_eigen_storage_structs`
            `num_rescalings_multipliers`
            `resource_index`
            `resource_preferences_flag`
            `resource_requirements_flag`
        """
        self._handle = None
        self._incarnated = False
        self._pattern_weight_list = None
        self._model_list = None
        self._num_leaves = None
        self._num_patterns = None
        self._num_states = None
        self._num_state_code_arrays = None
        self._num_partials = None
        self._num_model_matrices = None
        self._asrv_list = None
        self._num_prob_matrices = None
        self._num_eigen_storage_structs = None
        self._num_rescalings_multipliers = None

        self._pattern_weight_list = kwargs.get("pattern_weight_list")
        self._model_list = kwargs.get("model_list")
        self._num_leaves = kwargs.get("num_leaves")
        self.num_patterns = kwargs.get("num_patterns")
        self.num_states = kwargs.get("num_states")
        self._num_state_code_arrays = kwargs.get("num_state_code_arrays", 0)
        self._num_partials = kwargs.get("num_partials", 0)
        self.num_model_matrices = kwargs.get("num_model_matrices")
        self._asrv_list = kwargs.get("asrv_list")
        self._num_prob_matrices = kwargs.get("num_prob_matrices", 1)
        self._num_eigen_storage_structs = kwargs.get("num_eigen_storage_structs", 1)
        self._num_rescalings_multipliers = kwargs.get("num_rescalings_multipliers", 0)
        self._resource_index = kwargs.get("resource_index")
        self._resource_preferences_flag = kwargs.get("resource_preferences_flag", 0)
        self._resource_requirements_flag = kwargs.get("resource_requirements_flag", 0)

    def get_model_list(self):
        if self._model_list is not None:
            return tuple(self._model_list)
        nm = self.num_model_matrices
        return tuple([None]*nm)

    def set_model_list(self, v):
        if v is None:
            self._model_list = v
        else:
            if self._asrv_list is not None:
                if len(self._asrv_list) != len(self._model_list):
                    raise ValueError("If asrv_list and model_list are both specified then they must be the same length")
                for a, m in izip(self._asrv_list, v):
                    m.asrv = a
            self._num_model_matrices = len(v)
            self._model_list = tuple(v)
    model_list = property(get_model_list, set_model_list)

    def get_num_eigen_storage_structs(self):
        return self._num_eigen_storage_structs
    def set_num_eigen_storage_structs(self, v):
        if self._incarnated:
            raise RuntimeError("configuration attributes cannot be altered while the instance is incarnated!")
        i = int(v)
        if i < 0:
            raise ValueError("num_eigen_storage_structs cannot be less than 0")
        self._num_eigen_storage_structs = v
    num_eigen_storage_structs = property(get_num_eigen_storage_structs, set_num_eigen_storage_structs)

    def get_num_rescalings_multipliers(self):
        return self._num_rescalings_multipliers
    def set_num_rescalings_multipliers(self, v):
        if self._incarnated:
            raise RuntimeError("configuration attributes cannot be altered while the instance is incarnated!")
        i = int(v)
        if i < 0:
            raise ValueError("num_rescalings_multipliers cannot be less than 0")
        self._num_rescalings_multipliers = v
    num_rescalings_multipliers = property(get_num_rescalings_multipliers, set_num_rescalings_multipliers)

    def get_num_prob_matrices(self):
        return self._num_prob_matrices
    def set_num_prob_matrices(self, v):
        if self._incarnated:
            raise RuntimeError("configuration attributes cannot be altered while the instance is incarnated!")
        i = int(v)
        if i < 0:
            raise ValueError("num_prob_matrices cannot be less than 0")
        self._num_prob_matrices = v
    num_prob_matrices = property(get_num_prob_matrices, set_num_prob_matrices)

    def get_num_patterns(self):
        np = self._num_patterns
        if (np is None) and (self._pattern_weight_list is not None):
            return len(self._pattern_weight_list)
        return np
    def set_num_patterns(self, v):
        if self._incarnated:
            raise RuntimeError("configuration attributes cannot be altered while the instance is incarnated!")
        if v is None:
            self._num_patterns = v
            return
        i = int(v)
        if i < 0:
            raise ValueError("num_patterns cannot be less than 0")
        pw = self._pattern_weight_list
        if (pw is not None) and (len(pw) != i):
            raise ValueError("num_patterns must agree with the length of the list of pattern weights")
        self._num_patterns = v
    num_patterns = property(get_num_patterns, set_num_patterns)

    def get_pattern_weight_list (self):
        if self._pattern_weight_list is None:
            return ()
        return self._pattern_weight_list
    def set_pattern_weight_list(self, v):
        if self._incarnated:
            raise RuntimeError("configuration attributes cannot be altered while the instance is incarnated!")
        if v is None:
            self._pattern_weight_list = None
        else:
            self._pattern_weight_list = tuple(v)
            self._num_patterns = len(self._pattern_weight_list)
    pattern_weight_list = property(get_pattern_weight_list, set_pattern_weight_list)

    def get_num_partials(self):
        return self._num_partials
    def set_num_partials(self, v):
        if self._incarnated:
            raise RuntimeError("configuration attributes cannot be altered while the instance is incarnated!")
        i = int(v)
        if i < 0:
            raise ValueError("num_partials cannot be less than 0")
        self._num_partials = v
    num_partials = property(get_num_partials, set_num_partials)

    def get_resource_index(self):
        return self._resource_index
    def set_resource_index(self, v):
        if self._incarnated:
            raise RuntimeError("configuration attributes cannot be altered while the instance is incarnated!")
        i = int(v)
        if i < -1:
            raise ValueError("num_partials cannot be less than -1")
        self._resource_index = i
    resource_index = property(get_resource_index, set_resource_index)


    def get_num_state_code_arrays(self):
        return self._num_state_code_arrays
    def set_num_state_code_arrays(self, v):
        if self._incarnated:
            raise RuntimeError("configuration attributes cannot be altered while the instance is incarnated!")
        i = int(v)
        if i < 0:
            raise ValueError("num_state_code_arrays cannot be less than 0")
        self._num_state_code_arrays = v
    num_state_code_arrays = property(get_num_state_code_arrays, set_num_state_code_arrays)

    def get_num_leaves(self):
        return self._num_leaves
    def set_num_leaves(self, v):
        if self._incarnated:
            raise RuntimeError("configuration attributes cannot be altered while the instance is incarnated!")
        i = int(v)
        if i < 0:
            raise ValueError("num_leaves cannot be less than 0")
        self._num_leaves = v
    num_leaves = property(get_num_leaves, set_num_leaves)

    def get_num_states(self):
        return self._num_states
    def set_num_states(self, v):
        if self._incarnated:
            raise RuntimeError("configuration attributes cannot be altered while the instance is incarnated!")
        if v is None:
            self._num_states = v
            return
        i = int(v)
        if i < 2:
            raise ValueError("num_states cannot be less than 2")
        if self._model_list:
            for m in self._model_list:
                ns = m.num_states
                if (ns is not None) and ns != i:
                    raise ValueError("num_states must agree with num_states in all of the models is model_list") # could be relaxed by using multiple beagle instances
        self._num_states = v
    num_states = property(get_num_states, set_num_states)

    def get_asrv_list(self):
        if bool(self._asrv_list):
            return tuple(self._asrv_list)
        nm = self.num_model_matrices
        if nm is None:
            return tuple()
        if self._model_list:
            return tuple(i.asrv for i in self._model_list)
        return tuple()

    def set_asrv_list(self, v):
        if self._incarnated:
            raise RuntimeError("configuration attributes cannot be altered while the instance is incarnated!")
        if v is None:
            self._asrv_list = None
            return
        if self._model_list:
            if len(v) != len(self._model_list):
                raise ValueError("len(asrv_list) must equal len(model_list)") # could be relaxed by using multiple beagle instances
            for a, m in izip(v, self._model_list):
                m.asrv = a
        self._asrv_list = v
    asrv_list = property(get_asrv_list, set_asrv_list)

    def get_num_model_matrices(self):
        if self._num_model_matrices is None:
            if self._model_list is None:
                return None
            return len(self._model_list)
        return self._num_model_matrices
    def set_num_model_matrices(self, v):
        if self._incarnated:
            raise RuntimeError("configuration attributes cannot be altered while the instance is incarnated!")
        if v is None:
            self._num_model_matrices = None
        else:
            i = int(v)
            if i < 0:
                raise ValueError("num_model_matrices cannot be negative")
            if (self._model_list is None) or (i == len(self._model_list)):
                self._num_model_matrices = i
            else:
                raise ValueError("num_model_matrices and the length of the model list must agree (if both are used)")
    num_model_matrices = property(get_num_model_matrices, set_num_model_matrices)
    

    def _do_beagle_init(self):
        if self._incarnated:
            raise ValueError("Calculation instance has already been initialized. Duplicate intialization is not allowed")
        _LOG.debug("Calling cpytbeaglehon_init")
        asrv = self.asrv_list
        if self._resource_index is None:
            self._resource_index = 0
        if self._resource_preferences_flag is None:
            self._resource_preferences_flag = 0
        if self._resource_requirements_flag is None:
            self._resource_requirements_flag = 0

        if (self._model_list is not None) and len(self._model_list) > 0: 
            self.num_model_matrices = len(self._model_list)
            model_list = self._model_list
            self._model_list = list(model_list) # converting model_list from a tuple to a list
            models_wrappers_supplied = True
            self.num_states = model_list[0].num_states
        elif self.num_model_matrices is None:
            raise ValueError("A model list or num_model_matrices must be specified before initializing beagle")
        else:
            models_wrappers_supplied = False
            model_list = []

        if self.num_patterns is None:
            raise ValueError("A num_patterns list or pattern_weight_list must be specified before initializing beagle")
        resourceIndex = self._resource_index
        if resourceIndex is None:
            resourceIndex = -1
        arg_list = [self._num_leaves, 
                    self.num_patterns,
                    self.pattern_weight_list,
                    self.num_states,
                    self.num_state_code_arrays,
                    self.num_partials,
                    self.num_model_matrices,
                    asrv,
                    self.num_prob_matrices,
                    self.num_eigen_storage_structs,
                    self.num_rescalings_multipliers,
                    resourceIndex,
                    self._resource_preferences_flag,
                    self._resource_requirements_flag]
        _LOG.debug("Calling cpytbeaglehon_init with %s" %str(arg_list))
        self._handle = cpytbeaglehon_init(*arg_list)
        raw_models = cget_model_list(self._handle)
        for n, cmodel in enumerate(raw_models):
            try:
                a = asrv[n]
            except:
                a = None
            if models_wrappers_supplied:
                model_list[n]._reassign_environ(self, n, cmodel=cmodel, asrv=a)
            else:
                wrapped = DiscStateContTimeModel(cmodel=cmodel, num_states=self._num_states, model_index=n, calc_env=self, asrv=a)
                model_list.append(wrapped)
                
        if not models_wrappers_supplied:
            self._model_list = tuple(model_list)
        
        self._incarnated = True
        
        self._wrap_eigen_soln_structs = [EigenSolutionWrapper(index=n, like_calc_env=self) for n in range(self.num_eigen_storage_structs)]
        self._wrap_prob_mat = [ProbMatWrapper(index=n, like_calc_env=self) for n in range(self.num_prob_matrices)]
        self._wrap_partial = [PartialLikeWrapper(index=n, like_calc_env=self) for n in range(self.num_partials)]
        self._wrap_state_code_array = [StateCodeArrayWrapper(index=n, like_calc_env=self) for n in range(self.num_state_code_arrays)]
        self._wrap_rescalers = [RescalingMultiplier(index=n, like_calc_env=self) for n in range(self.num_rescalings_multipliers)]
        
        # caching is implemented by keeping a free, saved and calculated pool.
        #   if a new slot is needed, the free pool will be used, if it is empty
        #   then the calculated pool will be used. If that is also empty, then 
        #   a ValueError will be raised.
        self._eigen_cache = CalculatedCache(self._wrap_eigen_soln_structs, "EigenSolution")
        self._prob_mat_cache = CalculatedCache(self._wrap_prob_mat, "ProbMat")
        self._partial_cache = CalculatedCache(self._wrap_partial, "PartialLikelihood")
        self._state_code_cache = CalculatedCache(self._wrap_state_code_array, "StateCodeArray")
        self._rescalers_cache = CalculatedCache(self._wrap_rescalers, "Rescaler")

    def __del__(self):
        self.release_resources()

    def release_resources(self):
        if self._incarnated:
            _LOG.debug("Calling cpytbeaglehon_free")
            cpytbeaglehon_free(self._handle)
            self._model_list = ()
            self._handle = None
            self._incarnated = False
            
            del self._wrap_eigen_soln_structs
            del self._wrap_prob_mat
            del self._wrap_partial
            del self._wrap_state_code_array
            del self._wrap_rescalers
            del self._eigen_cache
            del self._prob_mat_cache
            del self._partial_cache
            del self._state_code_cache
            del self._rescalers_cache
            
            del self._free_eigen_storage_structs
            del self._saved_eigen_storage_structs
            del self._calculated_eigen_storage_structs
            del self._cached_eigen
    
            del self._free_prob_matrices
            del self._saved_prob_matrices
            del self._calculated_prob_matrices
            del self._cached_prob_matrices

    def calc_eigen_soln(self, model, model_state_hash, eigen_soln_caching=(CachingFacets.DO_NOT_SAVE,)):
        if not self._incarnated:
            self._do_beagle_init()
        _LOG.debug("LikeCalcEnvironment.calc_eigen_soln model=%d model_state_hash=%s eigen_soln_caching=%s" % (id(model), str(model_state_hash), str(eigen_soln_caching)))
        cf = eigen_soln_caching[0]
        e_cache = self._eigen_cache
        e_hash = EigenSolutionWrapper.calc_hash(model_state_hash)
        if cf == CachingFacets.RELEASE_THEN_SAVE:
            es_wrap = self._wrap_eigen_soln_structs[ eigen_soln_caching[1] ]
            if es_wrap.state_hash == e_hash:
                return es_wrap
            es_wrap = e_cache.make_writable(es_wrap)
        else:
            es_wrap = e_cache.get_from_cache(e_hash)

            if es_wrap is not None:
                if cf ==  CachingFacets.SAVE_ANYWHERE:
                    e_cache.incr_ref_count(es_wrap)
                return es_wrap
    
            es_wrap = e_cache.get_writable_object()

        try:
            es_wrap.calculate(model, model_state_hash)
        except:
            e_cache.release(es_wrap)
            raise
        
        if cf == CachingFacets.DO_NOT_SAVE:
            e_cache.flag_as_calculated(es_wrap)
        else:
            e_cache.save_obj(es_wrap)
        return es_wrap


    def calc_prob_from_eigen(self, edge_len, asrv, eigen_soln, prob_mat_caching=(CachingFacets.DO_NOT_SAVE,)):
        '''Returns a list of ProbMatWrapper objects."""
        '''
        if not self._incarnated:
            self._do_beagle_init()
        _LOG.debug("LikeCalcEnvironment.calc_prob_from_eigen edge_len=%f asrv=%s eigen_soln_index=%s eigen_state_id=%s prob_mat_caching=%s" % (float(edge_len), str(asrv), str(eigen_soln.index), str(eigen_soln.state_hash), str(prob_mat_caching)))
        cf = prob_mat_caching[0]

        if asrv is None:
            nc = 1
            rates = (1.0,)
            asrv_hash = _NONE_HASH
        else:
            nc = asrv.num_categories
            rates = asrv.rates
            asrv_hash = asrv.state_hash
        p_cache = self._prob_mat_cache
        to_return = []
        eigen_state_hash = eigen_soln.state_hash
        eff_edge_len_list = []
        eff_edge_len_hash_list = []
        pr_wrap_list = []
        try:
            for asrv_categ, rate in enumerate(rates):
                eff_edge_len = rate*edge_len
                edge_len_hash = repr(eff_edge_len)
                eff_edge_len_list.append(eff_edge_len)
                eff_edge_len_hash_list.append(edge_len_hash)
                curr_hash = ProbMatWrapper.calc_hash(eigen_state_hash, 
                                                     asrv_hash,
                                                     asrv_categ,
                                                     edge_len_hash)
                pr_wrap = p_cache.get_from_cache(curr_hash)
    
                if pr_wrap is not None:
                    # cache-hit...
                    if (cf == CachingFacets.RELEASE_THEN_SAVE) and (curr_hash != prob_mat_caching[1]):
                        raise ValueError("Calculation of probability matrix: RELEASE_THEN_SAVE specified, but the cache returned an unexpected object")
                else:
                    pr_wrap = p_cache.get_writable_object()
                pr_wrap_list.append(pr_wrap)

            ProbMatWrapper.calculate_list(pr_wrap_list,
                                          eigen_soln,
                                          asrv,
                                          eff_edge_len_list, 
                                          eigen_hash=eigen_state_hash,
                                          asrv_hash=asrv_hash, 
                                          eff_edge_len_hash=eff_edge_len_hash_list)
        except:
            for p in pr_wrap_list:
                p_cache.release(p)
            raise
        if cf == CachingFacets.DO_NOT_SAVE:
            for p in pr_wrap_list:
                p_cache.flag_as_calculated(p)
        else:
            for p in pr_wrap_list:
                p_cache.save_obj(p)
        return pr_wrap_list


    def get_prob_matrices(self, pr_wrap_list):
        _LOG.debug("LikeCalcEnvironment.get_prob_matrices([%s])" % (', '.join(['%d' % i.index for i in pr_wrap_list])))
        index_list = []
        for p in pr_wrap_list:
            assert(p.like_calc_env is self)
            assert(p.is_calculated)
            index_list.append(p.index)
        _LOG.debug("Calling cdsctm_get_pr_mats(%d, %s)" % (self._handle, str(index_list)))
        return cdsctm_get_pr_mats(self._handle, index_list)

    def set_state_code_array(self, leaf_index, leaf_data):
        if not self._incarnated:
            self._do_beagle_init()
        if not isinstance(leaf_data, tuple):
            leaf_data = tuple(leaf_data)
        o = self._wrap_state_code_array[leaf_index]
        self._state_code_cache.get_writable_object(o)
        try:
            _LOG.debug("Calling cdsctm_set_state_code(%s, leaf_index=%s, leaf_data=%s,..))" % (str(self._handle), str(leaf_index), str(leaf_data)))
            cdsctm_set_state_code(self._handle, leaf_index, leaf_data)
        except:
            self._state_code_cache.release(o)
            raise
        self._state_code_cache.save_obj(o)

    def tree_scorer(self, tree):
        if not self._incarnated:
            self._do_beagle_init()
        from pytbeaglehon.tree_scorer import TreeScorer
        return TreeScorer(like_calc_env=self, tree=tree)

    def start_partial_calculations(self, model):
        self._freeable_partials = {}
        self._freeable_prob_mats = {}
        self._queue_nd_to_prob_mats = {}
        self._num_queued_prob_mats = 0
        self._queue_nd_to_partial = {}
        self._num_queued_partials = 0
        self._in_partial_calcs = True
        self._num_prob_mats_avail_for_current = self.num_prob_matrices - len(self._saved_prob_matrices)
        self._queue_nd_order = []
        esi = model.eigen_soln_index # this triggers calculation of eigen solution

    def _force_calc_queued_partials(self, model):
        if self._num_queued_partials == 0:
            return
        assert(self._num_queued_partials == len(self._queue_nd_order))
        
        #TODO: this should be reorganized to reduces the number of Python->C->beagle calls!
        # for nd in self._queue_nd_order:
        #   do dependency pr mat calcs for each node
        # figure out what partial slots will be used for all calcs
        # for nd in self._queue_nd_order:
        #    add calc codes to list
        # call beagle calc partials
        partial_update_actions
        for nd in self._queue_nd_order:
            actions = self._queue_nd_to_partial[nd]
            assert(len(actions) == 3)
            left_child_action, right_child_action, curr_nd_action = actions
            left_child_dep, right_child_dep =  left_child_action[1][0], right_child_action[1][0]
            left_child_post, right_child_post =  left_child_action[1][1], right_child_action[1][1]
            if left_child_dep[0] is not None:
                assert(left_child_dep[0] == LikeCalcEnvironment._CALC_ACTION)
                left_child = left_child_dep[1]
                left_pmi = self.calc_prob_from_eigen(left_child.edge_length, 
                                                     asrv=model.asrv,
                                                     eigen_soln_index=model._eigen_soln_index,
                                                     eigen_state_id=model.state_hash,
                                                     prob_mat_caching=(SAVE_ANYWHERE,))
                left_child._LCE_prmat_index = left_pmi[0]
                left_child._LCE_edge_state_cache = left_pmi[1]
            if right_child_dep[0] is not None:
                assert(right_child_dep[0] == LikeCalcEnvironment._CALC_ACTION)
                right_child = right_child_dep[1]
                right_pmi = self.calc_prob_from_eigen(right_child.edge_length, 
                                                     asrv=model.asrv,
                                                     eigen_soln_index=model._eigen_soln_index,
                                                     eigen_state_id=model.state_hash,
                                                     prob_mat_caching=(SAVE_ANYWHERE,))
                right_child._LCE_prmat_index = right_pmi[0]
                right_child._LCE_edge_state_cache = right_pmi[1]

            assert(curr_nd_action[0] == LikeCalcEnvironment._CALC_ACTION)
            nd._LCE_is_internal = True
            if nd._LCE_save_partials:
                if nd._LCE_buffer_index is not None:
                    caching_arg = (CachingFacets.RELEASE_THEN_SAVE, nd._LCE_buffer_index)
                else:
                    caching_arg = (CachingFacets.SAVE_ANYWHERE)
            p = self._calc_single_partial(left_child, right_child)
            nd._LCE_buffer_index, nd._LCE_buffer_state_id = p
            
            if left_child_post[0] == LikeCalcEnvironment._FREE_ACTION:
                self._free_cached_prob_mat(left_pmi)
            if right_child_post[0] == LikeCalcEnvironment._FREE_ACTION:
                self._free_cached_prob_mat(right_pmi)
            if (left_child._LCE_is_internal) and (not left_child._LCE_save_partials):
                self._free_cached_partial(left_child._LCE_buffer_index)
                left_child._LCE_buffer_index = None
            if (right_child._LCE_is_internal) and (not right_child._LCE_save_partials):
                self._free_cached_partial(right_child._LCE_buffer_index)
                right_child._LCE_buffer_index = None
            
    def fetch_from_prob_mat_cache(self, model, node):
        assert(False)
        
        
    def _add_child_partial_calc(self, model, destination_node, child):
        pmi_state = self.fetch_from_prob_mat_cache(model, child)
        if pmi_state is None:
            if 1 + self._num_queued_prob_mats > self._num_prob_mats_avail_for_current:
                self._force_calc_queued_partials(model)
            self._num_queued_prob_mats += 1
            for_parent = self._queue_nd_to_prob_mats.setdefault(destination_node, [])
            t = ((LikeCalcEnvironment._CALC_ACTION, child), [LikeCalcEnvironment._FREE_ACTION, None])
            for_parent.append(t)
        else:
            ind, state_id = pmi_state
            if ind in self._saved_prob_matrices:
                self._saved_prob_matrices[ind][0] += 1
            else:
                self._saved_prob_matrices[ind] = [1, self._calculated_prob_matrices[ind]]
                del self._calculated_prob_matrices[ind]
            t = ((None,), [LikeCalcEnvironment._FREE_ACTION, pmi_state])
        return t

    def _partial_cached(self, nd, lse, rse):
        if lse is None or rse is None:
            return None
        sc = "%s %s" % (lse, rse)
        
        if sc == nd._LCE_partial_state_cache:
            try:
                i = self._cached_prob_matrices[sc]
                return (i, sc)
            except:
                return None
        return None
        
    def _add_binary_internal_node_to_partial_calc(self, model, destination_node, left_child, right_child):
        partial_cach = self._partial_cached(destination_node, left_child._LCE_edge_state_cache, right_child._LCE_edge_state_cache)
        if partial_cach is None:
            left_pmi = self._add_child_partial_calc(model, destination_node, left_child)
            right_pmi = self._add_child_partial_calc(model, destination_node, right_child)
            if 1 + self._num_queued_partials > self._num_partials_avail_for_current:
                self._force_calc_queued_partials(model)
            self._num_queued_partials += 1
            u = self._queue_nd_to_partial.setdefault(destination_node, [])
            u.append((LikeCalcEnvironment._DEPENDENCY_ACTION, left_pmi))
            u.append((LikeCalcEnvironment._DEPENDENCY_ACTION, right_pmi))
            u.append([LikeCalcEnvironment._CALC_ACTION, left_child, right_child])
            self._queue_nd_order.append(destination_node)
        else:
            ind, state_id = partial_cach
            if ind in self._saved_partials:
                self._saved_partials[ind][0] += 1
            else:
                self._saved_partials[ind] = [1, self._calculated_partials[ind]]
                del self._calculated_partials[ind]
        return pmi_state
        
    def add_internal_node_to_partial_calc(self, model, node):
        if not self._in_partial_calcs:
            raise ValueError("start_partial_calculations must be called before add_internal_node_to_partial_calc")
        c = node.children
        nc = len(c)
        if nc == 2:
            self._add_binary_internal_node_to_partial_calc(model, node, left_child=c[0], right_child=c[1])
        elif nc > 2:
            raise ValueError("Scoring trees with polytomies is not supported (yet).")
        elif nc == 1:
            raise ValueError("Scoring trees with nodes with outdegree = 1 is not supported (yet)")
        elif nc == 0:
            raise ValueError("add_internal_node_to_partial_calc should only be called for internal nodes")
            
    def end_partial_calculations(self, model):
        self._force_calc_queued_partials(model)
        if not self._in_partial_calcs:
            raise ValueError("start_partial_calculations must be called before add_internal_node_to_partial_calc")
        self._in_partial_calcs = False



def combine_state_id(*valist):
    return ' '.join([str(i) for i in valist])
##############################################################################
##  pytbeaglehon phylogenetic likelihood caluclations using beaglelib.
##
##  Copyright 2010 Mark T. Holder
##  All rights reserved.
##
##  See "LICENSE.txt" for terms and conditions of usage.
##
##############################################################################
