#! /usr/bin/env python
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
        if self.resource_index is None:
            raise ValueError("resource_index (the index of the computational resource) must be set before comp_resource_info can be accessed")
        return LikeCalcEnvironment.query_comp_resource_info(self.resource_index)
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
        self._initialized = False
        self.num_leaves = kwargs.get("num_leaves")
        self.num_patterns = kwargs.get("num_patterns")
        self.pattern_weight_list = kwargs.get("pattern_weight_list", [])
        self.num_states = kwargs.get("num_states")
        self.num_state_code_arrays = kwargs.get("num_state_code_arrays")
        self.num_partials = kwargs.get("num_partials")
        self.num_inst_rate_matrices = kwargs.get("num_inst_rate_matrices")
        self.asrv_list = kwargs.get("asrv_list", [])
        self.num_prob_matrices = kwargs.get("num_prob_matrices")
        self.num_eigen_storage_structs = kwargs.get("num_eigen_storage_structs")
        self.num_rescalings_multipliers = kwargs.get("num_rescalings_multipliers")
        self.resource_index = kwargs.get("resource_index")
        self.resource_preferences_flag = kwargs.get("resource_preferences_flag")
        self.resource_requirements_flag = kwargs.get("resource_requirements_flag")
        self._char_edge_model = ()
    def _do_beagle_init(self):
        if self._initialized:
                raise ValueError("Calculation instance has already been initialized. Duplicate intialization is not allowed")
        _LOG.debug("Calling cpytbeaglehon_init")
        if self.resource_index is None:
            self.resource_index = 0
        if self.resource_preferences_flag is None:
            self.resource_preferences_flag = 0
        if self.resource_requirements_flag is None:
            self.resource_requirements_flag = 0
        self._handle = cpytbeaglehon_init( self.num_leaves, 
                            self.num_patterns,
                            self.pattern_weight_list,
                            self.num_states,
                            self.num_state_code_arrays,
                            self.num_partials,
                            self.num_inst_rate_matrices,
                            self.asrv_list,
                            self.num_prob_matrices,
                            self.num_eigen_storage_structs,
                            self.num_rescalings_multipliers,
                            self.resource_index,
                            self.resource_preferences_flag,
                            self.resource_requirements_flag)
        raw_models = cget_model_list(self._handle)
        self._char_edge_model = ()
        w = []
        for n, i in enumerate(raw_models):
            try:
                a = self.asrv_list[n]
            except:
                a = None
            wrapped = DiscStateContTimeModel(cmodel=i, num_states=self.num_states, model_index=n, calc_env=self, asrv=a)
            w.append(wrapped)
        self._char_edge_model = tuple(w)
        
        self._initialized = True

    def __del__(self):
        self.release_resources()

    def release_resources(self):
        if self._initialized:
            _LOG.debug("Calling cpytbeaglehon_free")
            cpytbeaglehon_free(self._handle)
            self._char_edge_model = ()
            self._handle = None
            self._initialized = False
