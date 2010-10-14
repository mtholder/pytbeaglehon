/**
  pytbeaglehon phylogenetic likelihood caluclations using beaglelib.

  Copyright 2010 Mark T. Holder (except where noted). All rights reserved.

  See "LICENSE.txt" for terms and conditions of usage.
  
  
  
  Code for writing extensions was adapted from Alex Martelli's example code 
    in the Python Cookbook 17.1

*/
#include "py_beagle.h"
#include "py_asrv.h"
#include "calc_instance.h"


/* exception class (will be assigned in initdsct_model) */
static PyObject * CLAUnderflowError;




static PyObject* cpytbeaglehon_init(PyObject *self, PyObject *args) {
	int numLeaves;
	long numPatterns;
	int num_states;
	int num_rate_categories;
	int num_state_code_arrays;
	int numPartialStructs;
	int num_model_matrices;
	int numProbMats;
	int numEigenStorage;
	int num_rescalings_multipliers;
	int resource_arg;
	long resource_flag;
	long handle;
	if (!PyArg_ParseTuple(args, "iliiiiiiiiil", &numLeaves,
											  &numPatterns,
											  &num_states,
											  &num_rate_categories,
											  &num_state_code_arrays,
											  &numPartialStructs,
											  &num_model_matrices,
											  &numProbMats,
											  &numEigenStorage,
											  &num_rescalings_multipliers,
											  &resource_arg,
											  &resource_flag)) {
		return 0L;
	}
	if (numLeaves < 0) {
		PyErr_SetString(PyExc_ValueError, "The number of leaves cannot be negative");
		return 0;
	}
	if (numPatterns < 0) {
		PyErr_SetString(PyExc_ValueError, "The number of patterns cannot be negative");
		return 0;
	}
	if (num_states < 2) {
		PyErr_SetString(PyExc_ValueError, "The number of states cannot be less than 1");
		return 0;
	}
	if (num_rate_categories < 1) {
		PyErr_SetString(PyExc_ValueError, "The number of rate categories cannot be less than 1");
		return 0;
	}
	if (num_state_code_arrays < 0) {
		PyErr_SetString(PyExc_ValueError, "The number of state code arrays cannot be negative");
		return 0;
	}
	if (numPartialStructs < 0) {
		PyErr_SetString(PyExc_ValueError, "The number of partial likelihood arrays cannot be negative");
		return 0;
	}
	if (num_model_matrices < 0) {
		PyErr_SetString(PyExc_ValueError, "The number of model matrices cannot be negative");
		return 0;
	}
	if (numProbMats < 0) {
		PyErr_SetString(PyExc_ValueError, "The number of transition probability matrices cannot be negative");
		return 0;
	}
	if (numEigenStorage < 0) {
		PyErr_SetString(PyExc_ValueError, "The number of eigen solutions cannot be negative");
		return 0;
	}
	if (num_rescalings_multipliers < 0) {
		PyErr_SetString(PyExc_ValueError, "The number of rescaling buffers cannot be negative");
		return 0;
	}
	handle = create_likelihood_calc_instance(numLeaves,
											 numPatterns,
											 num_states,
											 num_rate_categories,
											 num_state_code_arrays,
											 numPartialStructs,
											 num_model_matrices,
											 numProbMats,
											 numEigenStorage,
											 num_rescalings_multipliers,
											 resource_arg,
											 resource_flag);
	if (handle < 0) {
		PYTBEAGLEHON_DEBUG_PRINTF("Could not allocate likelihood_calculator_instance\n");
		PyErr_NoMemory();
		return 0L;
	}
	return PyInt_FromLong(handle);
}




static PyMethodDef dsct_model_module_functions[] = {
#if 0
    {"get_full_la_vals", get_full_la_vals, METH_VARARGS,
		"Takes a FullLAObj returns a list of pattern lnL, a list of frequencies for categ x state, the \"root\" conditional likelihoods, and a list of rescalers (either ints or floats).\nCLA's are arranged site x categ x state"},
	{"get_cla_vals", get_cla_vals, METH_VARARGS,
		"Takes a cla object and returns a list of conditional likelihoods, and a list of rescalers (either ints or floats).\nCLA's are arranged site x categ x state"},
	{"calc_ln_L_across_term", ln_L_terminal_edge,	 METH_VARARGS,
		"Calculates ln  likelihood across a terminal branch: args (full LA obj, tip data, tip pmat, par cla, beg cat, end cat, rescale_threshold); Returns ln Likelihood."},
	{"calc_ln_L_across_internal", ln_L_internal_edge,	 METH_VARARGS,
		"Calculates ln  likelihood across an internal branch: args (full LA obj, child cla, pmat, par cla, beg cat, end cat, rescale_threshold); Returns ln Likelihood."},
	{"calc_anc_from_two_tips", two_tip_cla,	 METH_VARARGS,
		"Calculates conditional likelihood array from two tips: args (par cla obj, first tip data, first pmat, second tip data, second pmat, beg cat, end cat, rescale_threshold); Returns rescale_threshold."},
	{"calc_anc_from_one_tip_one_intern", one_tip_one_internal_cla,	 METH_VARARGS,
		"Calculates conditional likelihood array from 1 tip and one internal: args (par cla obj,  tip data, tip pmat, child cla, internal pmat, beg cat, end cat, rescale_threshold); Returns rescale_threshold."},
	{"calc_anc_from_two_internals", two_internal_cla,	 METH_VARARGS,
		"Calculates conditional likelihood array from two tips: args (par cla obj, first child cla, first pmat, second child cla, second pmat, beg cat, end cat, rescale_threshold); Returns rescale_threshold."},
	{"add_leaf_to_cla", add_leaf_to_cla,	 METH_VARARGS,
		"Conditions cla on another tip: args (par cla obj,  tip data,  pmat, beg cat, end cat, rescale_threshold); Returns rescale_threshold."},
	{"add_internal_to_cla", add_internal_to_cla,	 METH_VARARGS,
		"Conditions cla on another internal: args (par cla obj,  child cla,  pmat, beg cat, end cat, rescale_threshold); Returns rescale_threshold."},
	{"cleaf_data_ctor", cleaf_data_ctor, METH_VARARGS,
		"initializer -- takes a list of ints and StateSetLookupStruct"},
	{"cla_ctor", cla_ctor, METH_VARARGS,
		"initializer -- takes n_sites, nStates, n_categ"},
	{"full_la_set_freqs", full_la_set_freqs, METH_VARARGS,
		"Takes a full_la object, a list of category frequencies and a list of lists of state frequencies"},
	{"full_la_ctor", full_la_ctor, METH_VARARGS,
		"initializer -- takes n_sites, nStates, n_categ"},
	{"cdsctm_ctor", cdsctm_ctor,	 METH_VARARGS,
		"initializer -- takes the number of states (integer)."},
	{"cdsctm_get_dimen", cdsctm_get_dimen, METH_VARARGS,
		"returns the number of states"},
	{"cdsctm_get_q_mat", cdsctm_get_q_mat, METH_VARARGS,
		"Returns the Q-matrix as a list of lists of floats"},
	{"cdsctm_get_q_mat_row", cdsctm_get_q_mat_row, METH_VARARGS,
		"Returns a row the Q-matrix as a list of floats (takes row index)."},
	{"cdsctm_get_q_mat_cell", cdsctm_get_q_mat_cell, METH_VARARGS,
		"Returns a cell of the Q-matrix (takes row and col. indices)"},
	{"cdsctm_set_q_mat", cdsctm_set_q_mat, METH_VARARGS,
		"Replaces the Q-matrix (takes a list of list of floats)."},
	{"cdsctm_set_q_mat_row", cdsctm_set_q_mat_row, METH_VARARGS,
		"Replaces a row of the Q-matrix (takes a row index and list)"},
	{"cdsctm_set_q_mat_cell", cdsctm_set_q_mat_cell, METH_VARARGS,
		"Replaces a row of the Q-matrix (takes a row and col. index and value)"},


	{"get_num_computational_resources", get_num_computational_resources, METH_VARARGS,
		"Queries installed libraries for the number of possible likelihood computational configurations"},
	{"get_resource_info", get_resource_info, METH_VARARGS,
		"Takes an integer, and returns a tuple of (resourceNumber, flags, name, implName) for that resource"},
	{"cpytbeaglehon_init", cpytbeaglehon_init, METH_VARARGS,
		"Creates an new instance of a likelihood calculation context (allocates storage, etc).  Returns a handle that must be used to provide an \"address space\" for the indices in most other cPhyProb method calls."},
	{"cphyprob_free", cphyprob_free, METH_VARARGS,
		"Frees memory associated with a LikeCalculatorInstance"},


	/* The following methods REQUIRE a cphyprob calc context instance from cpytbeaglehon_init as the first argument */
	{"calc_pmat", ccalc_pmat,	 METH_VARARGS,
		"Calculates a P-Matrix from a Q-matrix and branch length.\nTakes a pmat_array_type, dsct_model_type, and branch length"},
	{"calc_pmat_array", ccalc_pmat_array,	 METH_VARARGS,
		"Calculates a P-Matrix from a Q-matrix and branch length.\nTakes a pmat_array_type, list of dsct_model_type, and list of branch length"},
	{"cpmat_array_ctor", cpmat_array_ctor,	 METH_VARARGS,
		"Initializer -- takes the # of matrices and # of states (positive integers)."},
	{"cpmat_array_get_n_matrices", cpmat_array_get_n_matrices, METH_VARARGS,
		"Returns the number of matrices in a PMatArrayObj"},
	{"cpmat_array_get_n_states", cpmat_array_get_n_states, METH_VARARGS,
		"Returns the number of states"},
	{"cpmat_array_get_mat_list", cpmat_array_get_mat_list, METH_VARARGS,
		"Returns a list P-matrix as a list of  list of lists of floats"},
	{"cpmat_array_get_mat", cpmat_array_get_mat, METH_VARARGS,
		"Takes an index and returns a P-matrix as a list of lists of floats"},


	/* The following methods do NOT require a cphyprob calc context instance */
	{"cstate_set_lookup_ctor", cstate_set_lookup_ctor,	 METH_VARARGS,
		"initializer -- takes a list of lists.  Returns a StateSetLookupStruct."},
#endif

	{"casrvo_ctor", casrvo_ctor, METH_VARARGS,
		"intializer -- takes shape parameter, ncat, dcst_model.RateHetType facet."},
	{"casrvo_get_n_cat", casrvo_get_n_cat, METH_VARARGS,
		"Returns the number of categories."},
	{"casrvo_get_rates", casrvo_get_rates, METH_VARARGS,
		"Returns a list of the rate multipliers."},
	{"casrvo_get_shape", casrvo_get_shape, METH_VARARGS,
		"Returns the shape parameter."},
	{"casrvo_set_shape", casrvo_set_shape, METH_VARARGS,
		"Sets the shape parameter (should be a positive double)."},

    {0, 0}
};





/* module entry-point (module-initialization) function */
void initdisc_state_cont_time_model(void)
{
	PyObject * d, * m;
		/* Create the module, with its functions */
	m = Py_InitModule("disc_state_cont_time_model", dsct_model_module_functions);

	asrv_type.ob_type = &PyType_Type;

	d = PyModule_GetDict(m);

	CLAUnderflowError = PyErr_NewException("claUnderflow.error", NULL, NULL);
	PyDict_SetItemString(d, "error", CLAUnderflowError);

}

