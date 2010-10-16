#! /usr/bin/env python
import itertools
from pytbeaglehon.ccore.disc_state_cont_time_model import cpytbeaglehon_init, cpytbeaglehon_free, cget_num_comp_resources, cget_comp_resource_info, cget_model_list
from pytbeaglehon import DiscStateContTimeModel
from pytbeaglehon import get_logger

_LOG = get_logger(__name__)

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




class LikeCalcEnvironment(object):
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
                for a, m in itertools.izip(self._asrv_list, v):
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
            return tuple([i.asrv] for i in self._model_list)
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
            for a, m in itertools.izip(v, self._model_list):
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
            if (self._model_list is None) or (i != len(self._model_list)):
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
        self._handle = cpytbeaglehon_init( self._num_leaves, 
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
                            self._resource_requirements_flag)
        raw_models = cget_model_list(self._handle)
        self._model_list = ()
        for n, i in enumerate(raw_models):
            try:
                a = asrv[n]
            except:
                a = None
            if models_wrappers_supplied:
                self._model_list[n]._reassign_environ(self, n, cmodel=cmodel, asrv=a)
            else:
                wrapped = DiscStateContTimeModel(cmodel=i, num_states=self._num_states, model_index=n, calc_env=self, asrv=a)
                model_list.append(wrapped)
                
        if not models_wrappers_supplied:
            self._model_list = tuple(model_list)
        
        self._incarnated = True

    def __del__(self):
        self.release_resources()

    def release_resources(self):
        if self._incarnated:
            _LOG.debug("Calling cpytbeaglehon_free")
            cpytbeaglehon_free(self._handle)
            self._model_list = ()
            self._handle = None
            self._incarnated = False
