#include "py_discrete_state_model.h"
#include "discrete_state_model.h"
#include "py_util.h"


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
	dim = dsct_model_obj->dim;
	to_return = tupleToDoubleMatrix(tuple_obj, dsct_model_obj->qMat, dim, dim, 1);
	dsct_model_obj->eigenCalcIsDirty = 1;
	return to_return;
}




PyObject* cdsctm_calc_eigens(PyObject *self, PyObject *args) {
	int eigenIndex;
	PyObject * dsct_model_py_obj;
	DSCTModelObj *dsct_model_obj;
	if (!PyArg_ParseTuple(args, "O!O!", &dsct_model_type, &dsct_model_py_obj, &eigenIndex))
		return 0L;
	dsct_model_obj = (DSCTModelObj *)(dsct_model_py_obj);
	dsct_model_obj->eigenBufferIndex = eigenIndex;
	dsct_model_obj->eigenCalcIsDirty = 1;
	if (recalc_eigen_mat(dsct_model_obj) == 0)
	    return 0L;
	return none();
}
