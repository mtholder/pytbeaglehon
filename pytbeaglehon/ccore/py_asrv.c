#include "py_asrv.h"
#include "asrv.h"
#include "py_util.h"



/**
  pytbeaglehon phylogenetic likelihood caluclations using beaglelib.

  Copyright 2010 Mark T. Holder (except where noted). All rights reserved.

  See "LICENSE.txt" for terms and conditions of usage.
  
  
  
  Code for writing extensions was adapted from Alex Martelli's example code 
    in the Python Cookbook 17.1

*/



/* Forward declare the python type object*/
PyTypeObject asrv_type;

/*******************************************************************************
 * Python type wrappers
 */

PyTypeObject asrv_type = {
    PyObject_HEAD_INIT(0)	  	/* initialize to 0 to ensure Win32 portability  */
    0,						  	/* ob_size */
    "asrv_obj",	/* tp_name */
    sizeof(ASRVObj),		/* tp_basicsize */
    0,						  	/* tp_itemsize */
    /* methods */
    (destructor)asrv_obj_dtor, /* tp_dealloc */
    /* implied by ISO C: all zeros thereafter, i.e., no other method */
};

/*******************************************************************************
 * ASRV adaptor functions
 */

/** the factory function */
PyObject* casrvo_ctor(PyObject *self, PyObject *args) {
	int n_categ, mode_enum;
	double val;
	unsigned dim;
	if (!PyArg_ParseTuple(args, "dii", &val, &n_categ, &mode_enum))
		return 0;
	if (n_categ < 1) {
		PyErr_SetString(PyExc_ValueError, "The number of categories parameter must be greater than 0");
		return 0;
	}
	dim = (unsigned) n_categ;
	return (PyObject*) asrv_obj_new(dim, mode_enum, val);
}

PyObject* casrvo_get_shape(PyObject *self, PyObject *args) {
	PyObject *asrh_py_obj;
	ASRVObj * asrv;
	if (!PyArg_ParseTuple(args, "O!", &asrv_type, &asrh_py_obj))
		return 0L;
	asrv = (ASRVObj *)(asrh_py_obj);
	return PyFloat_FromDouble(asrv->param);
}
PyObject* casrvo_get_n_cat(PyObject *self, PyObject *args) {
	PyObject *asrh_py_obj;
	ASRVObj * asrv;
	if (!PyArg_ParseTuple(args, "O!", &asrv_type, &asrh_py_obj))
		return 0L;
	asrv = (ASRVObj *)(asrh_py_obj);
	return PyInt_FromLong((long)(asrv->n));
}
PyObject* casrvo_set_shape(PyObject *self, PyObject *args) {
	double val;
	PyObject *mfo_py_obj;
	if (!PyArg_ParseTuple(args, "O!d", &asrv_type, &mfo_py_obj, &val))
		return 0L;
	if (val <= 0.0) {
		PyErr_SetString(PyExc_ValueError, "The shape parameter must be > 0.0");
		return 0L;
	}
	internal_asrv_set_shape((ASRVObj *)(mfo_py_obj), val);
	return none();
}
PyObject* casrvo_get_rates(PyObject *self, PyObject *args) {
	PyObject *asrh_py_obj;
	ASRVObj * asrv;
	if (!PyArg_ParseTuple(args, "O!", &asrv_type, &asrh_py_obj))
		return 0L;
	asrv = (ASRVObj *)(asrh_py_obj);
	return doubleArrayToList(asrv->val, asrv->n);
}
