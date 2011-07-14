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
#   if defined(API_TRACE_PRINTING) && API_TRACE_PRINTING
        const struct LikeCalculatorInstance *LCI = getLikeCalculatorInstance(0);
    	if (dim == 4) {
    	    PYTBEAGLEHON_DEBUG_PRINTF2("/* cAPI Call */ const double row0[4] = {%.10lf, %.10lf, %.10lf, %.10lf};\n", dsct_model_obj->qMat[0][0], dsct_model_obj->qMat[0][1], dsct_model_obj->qMat[0][2], dsct_model_obj->qMat[0][3]);
    	    PYTBEAGLEHON_DEBUG_PRINTF2("/* cAPI Call */ const double row1[4] = {%.10lf, %.10lf, %.10lf, %.10lf};\n", dsct_model_obj->qMat[1][0], dsct_model_obj->qMat[1][1], dsct_model_obj->qMat[1][2], dsct_model_obj->qMat[1][3]);
    	    PYTBEAGLEHON_DEBUG_PRINTF2("/* cAPI Call */ const double row2[4] = {%.10lf, %.10lf, %.10lf, %.10lf};\n", dsct_model_obj->qMat[2][0], dsct_model_obj->qMat[2][1], dsct_model_obj->qMat[2][2], dsct_model_obj->qMat[2][3]);
    	    PYTBEAGLEHON_DEBUG_PRINTF2("/* cAPI Call */ const double row3[4] = {%.10lf, %.10lf, %.10lf, %.10lf};\n", dsct_model_obj->qMat[3][0], dsct_model_obj->qMat[3][1], dsct_model_obj->qMat[3][2], dsct_model_obj->qMat[3][3]);
    	    PYTBEAGLEHON_DEBUG_PRINTF2("/* cAPI Call */ const double * mat[4] = {&row0[0], &row1[0], &row2[0], &row3[0]};\n");
            PYTBEAGLEHON_DEBUG_PRINTF2("/* cAPI Call */ setQMatForHandle(%d, %d, &mat[0]);\n", LCI->handleForSelf getTraceModeModelIndex(LCI, dsct_model_obj));
        }
        else {
            PYTBEAGLEHON_DEBUG_PRINTF2("/* cAPI Call */ setQMatForHandle(%d, %d, QMAT_GOESHERE);\n", LCI->handleForSelf getTraceModeModelIndex(LCI, dsct_model_obj));
        }
#   endif
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
	dsct_model_obj = (DSCTModelObj *)(dsct_model_py_obj);
#   if defined(API_TRACE_PRINTING) && API_TRACE_PRINTING
        const struct LikeCalculatorInstance *LCI = getLikeCalculatorInstance(0);
        PYTBEAGLEHON_DEBUG_PRINTF1("/* cAPI Call */ dsct_model_obj = LCI->probModelArray[%d];\n", getTraceModeModelIndex(LCI, dsct_model_obj));
    	PYTBEAGLEHON_DEBUG_PRINTF1("/* cAPI Call */ dsct_model_obj->eigenBufferIndex=%d; dsct_model_obj->eigenCalcIsDirty = 1;\n", eigenIndex);
    	PYTBEAGLEHON_DEBUG_PRINTF("/* cAPI Call */ rc = recalc_eigen_mat(dsct_model_obj); if (rc == 0) {fprintf(stderr, \"recalc_eigen_mat failed\"); return 1;}\n");
#   endif
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
    const struct LikeCalculatorInstance * lci;
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
#   if defined(DEBUG_PRINTING) && DEBUG_PRINTING
        unsigned i;
#   endif
    long handle;
    int eigenIndex;
    unsigned numToCalc;
	PyObject * edge_len_list_obj, * pr_mat_ind_list_obj;
    const struct LikeCalculatorInstance * lci;
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
#   if defined(API_TRACE_PRINTING) && API_TRACE_PRINTING
	    PYTBEAGLEHON_DEBUG_PRINTF3("/* cAPI Call */ handle=%ld; eigenIndex=%d; numToCalc=%u;", handle, eigenIndex, numToCalc);
	    for (i = 0; i < numToCalc; ++i) {
    	    PYTBEAGLEHON_DEBUG_PRINTF2("LCI->edgeLenScratch[%d] = %.10lf; ", i, lci->edgeLenScratch[i]);
        }
	    PYTBEAGLEHON_DEBUG_PRINTF(";\n/* cAPI Call */ ");
	    for (i = 0; i < numToCalc; ++i) {
    	    PYTBEAGLEHON_DEBUG_PRINTF2("LCI->probMatIndexScratch[%d] = %d; ", i, lci->probMatIndexScratch[i]);
        }
	    PYTBEAGLEHON_DEBUG_PRINTF("\n/* cAPI Call */ rc = calcPrMatsForHandle(handle, eigenIndex, numToCalc, LCI->edgeLenScratch, LCI->probMatIndexScratch); if (rc != 0) {return rc;}\n");
#   endif
	if (calcPrMatsForHandle(handle, eigenIndex, numToCalc, lci->edgeLenScratch, lci->probMatIndexScratch) != BEAGLE_SUCCESS) {
	    PyErr_SetString(PyExc_RuntimeError, "calcPrMatsForHandle call failed");
	    return 0L;
	}
	return none();
}
