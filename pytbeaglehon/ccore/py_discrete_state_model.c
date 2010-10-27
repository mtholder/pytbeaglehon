#include <libhmsbeagle/beagle.h>
#include "py_discrete_state_model.h"
#include "discrete_state_model.h"
#include "py_util.h"
#include "internal_like_calc_env.h"

/* Q matrix setter */
PyObject* cdsctm_set_q_mat(PyObject *self, PyObject *args) {
	unsigned dim;
	PyObject * tuple_obj = 0L;
	PyObject * dsct_model_py_obj;
	PyObject * to_return;
	DSCTModelObj *dsct_model_obj;
	if (!PyArg_ParseTuple(args, "O!O!", &dsct_model_type, &dsct_model_py_obj, &PyTuple_Type, &tuple_obj))
		return 0L;
	dsct_model_obj = (DSCTModelObj *)(dsct_model_py_obj);
	PYTBEAGLEHON_DEBUG_PRINTF1("Setting qmat for model at %ld\n", (long) dsct_model_obj);
	dim = dsct_model_obj->dim;
	to_return = tupleToDoubleMatrix(tuple_obj, dsct_model_obj->qMat, dim, dim, 1);
	dsct_model_obj->eigenCalcIsDirty = 1;
	PYTBEAGLEHON_DEBUG_PRINTF1("qmat for model at %ld set.\n", (long) dsct_model_obj);
	return to_return;
}

PyObject* cdsctm_calc_eigens(PyObject *self, PyObject *args) {
	int eigenIndex;
	PyObject * dsct_model_py_obj;
	DSCTModelObj * dsct_model_obj;
	if (!PyArg_ParseTuple(args, "O!i", &dsct_model_type, &dsct_model_py_obj, &eigenIndex))
		return 0L;
	PYTBEAGLEHON_DEBUG_PRINTF2("In C. calling cdsctm_calc_eigens(%ld, %d)\n", (long) dsct_model_py_obj, eigenIndex);
	dsct_model_obj = (DSCTModelObj *)(dsct_model_py_obj);
	dsct_model_obj->eigenBufferIndex = eigenIndex;
	dsct_model_obj->eigenCalcIsDirty = 1;
	if (recalc_eigen_mat(dsct_model_obj) == 0)
	    return 0L;
	return none();
}


PyObject* cdsctm_get_pr_mats(PyObject *self, PyObject *args) {
    long handle;
    unsigned i, numToCalc;
	PyObject * pr_mat_ind_list_obj;
    struct LikeCalculatorInstance * lci;
    PyObject *lp;
	PyObject *el_obj;

	if (!PyArg_ParseTuple(args, "lO!", &handle, &PyList_Type, &pr_mat_ind_list_obj))
		return 0L;
    lci = getLikeCalculatorInstance(handle);
    if (lci == 0L) {
		PyErr_SetString(PyExc_IndexError, "LikeCalculatorInstance handle out of range");
		return 0L;
    }
	if (listToUnsignedArrayMaxSize(pr_mat_ind_list_obj, lci->probMatIndexScratch, lci->numProbMats, &numToCalc) == 0) {
		PyErr_SetString(PyExc_IndexError, "Number of probability matrix exceeds the number that exist!");
		return 0L;
    }    
    
	lp = PyList_New(numToCalc);
	if (lp == 0L)
		return 0L;
    
    for (i = 0; i < numToCalc; ++i) {
        if (fetchPrMat(handle, lci->probMatIndexScratch[i], lci->probMatScratch[0]) != BEAGLE_SUCCESS) {
            PyErr_SetString(PyExc_RuntimeError, "fetchPrMat call failed");
            goto errorExit;
        }
		el_obj = doubleMatToList((const double**)lci->probMatScratch, lci->numStates, lci->numStates);
		if (el_obj == 0L) {
			goto errorExit;
		}
		PyList_SetItem(lp, i, el_obj);
    }    
    return lp;
    errorExit:
        Py_DECREF(lp);
        return 0L;
}
PyObject* cdsctm_calc_pr_mats(PyObject *self, PyObject *args) {
    long handle;
    int eigenIndex;
    unsigned numToCalc;
	PyObject * edge_len_list_obj, * pr_mat_ind_list_obj;
    struct LikeCalculatorInstance * lci;
	PYTBEAGLEHON_DEBUG_PRINTF("Entering cdsctm_calc_pr_mats\n");
	if (!PyArg_ParseTuple(args, "liO!O!", &handle, &eigenIndex, &PyList_Type, &edge_len_list_obj, &PyList_Type, &pr_mat_ind_list_obj))
		return 0L;
	PYTBEAGLEHON_DEBUG_PRINTF("getting lci in cdsctm_calc_pr_mats\n");
    lci = getLikeCalculatorInstance(handle);
    if (lci == 0L) {
		PyErr_SetString(PyExc_IndexError, "LikeCalculatorInstance handle out of range");
		return 0L;
    }
	PYTBEAGLEHON_DEBUG_PRINTF("unpacking edge lenths in  cdsctm_calc_pr_mats\n");
	if (listToDoubleArrayMaxSize(edge_len_list_obj, lci->edgeLenScratch, lci->numProbMats, &numToCalc) == 0)
	    return 0L;
	PYTBEAGLEHON_DEBUG_PRINTF("unpacking prob mat buffers in  cdsctm_calc_pr_mats\n");
	if (listToUnsignedArray(pr_mat_ind_list_obj, lci->probMatIndexScratch, numToCalc) == 0) {
		PyErr_SetString(PyExc_IndexError, "edge length list and prob mat index list must be the same length");
	}
	PYTBEAGLEHON_DEBUG_PRINTF("Calling calcPrMats...");
	if (calcPrMats(handle, eigenIndex, numToCalc, lci->edgeLenScratch, lci->probMatIndexScratch) != BEAGLE_SUCCESS) {
	    PyErr_SetString(PyExc_RuntimeError, "calcPrMats call failed");
	    return 0L;
	}
	return none();
}
