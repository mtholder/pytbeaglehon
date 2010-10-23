#include "py_util.h"


/**
  pytbeaglehon phylogenetic likelihood caluclations using beaglelib.

  Copyright 2010 Mark T. Holder (except where noted). All rights reserved.

  See "LICENSE.txt" for terms and conditions of usage.
  
  
  
  Code for writing extensions was adapted from Alex Martelli's example code 
    in the Python Cookbook 17.1

*/

/*******************************************************************************
 * Function definitions
 */

		/***********************************************************************
		 * Utility functions for Python binding
		 */

PyObject* none(void) {
	Py_INCREF(Py_None);
	return Py_None;
}

/**
 * (internal) converts an python index on [-dim,dim) to the positive version
 * of the index.
 * Returns a non-negative index for success or -1.
 * if -1 is returned then the function will have generated a
 * IndexError("list index out of range") already.
 *
 */
int pyIndToCInd(int py_ind, unsigned dim) {
	PYTBEAGLEHON_DEBUG_PRINTF2("In pyIndToCInd(%d,%u)\n", py_ind, dim);
	if (py_ind < 0)
		py_ind += (int)dim;
	if (py_ind < (int)dim)
		return py_ind;
	PyErr_SetString(PyExc_IndexError, "list index out of range");
	return -1;
}
/**
 * (internal) Converts a double array to a python list of floats
 *
 *	\assert(arr != 0L || len == 0)
 */
PyObject * doubleArrayToList(const double *arr, unsigned len) {
	PyObject *lp;
	PyObject *el_obj;
	unsigned i;
	assert(arr != 0L || len == 0);
	lp = PyList_New(len);
	if (lp == 0L)
		return 0L;
	for (i = 0; i < len; ++i) {
		el_obj = PyFloat_FromDouble(arr[i]);
		if (el_obj == 0L) {
			Py_DECREF(lp);
			return 0L;
		}
		PyList_SetItem(lp, i, el_obj);
	}
	return lp;
}
/**
 * (internal)  Converts a rectangular 2D matrix of floats to a python list of
 * 	list of floats.
 *
 *	\assert(arr != 0L || len == 0)
 */
PyObject * doubleMatToList(const double **arr, unsigned n_rows, unsigned n_cols) {
	PyObject *lp;
	PyObject *el_obj;
	unsigned i;
	assert(arr != 0L || n_rows == 0);
	lp = PyList_New(n_rows);
	if (lp == 0L)
		return 0L;
	for (i = 0; i < n_rows; ++i) {
		el_obj = doubleArrayToList(arr[i], n_cols);
		if (el_obj == 0L) {
			Py_DECREF(lp);
			return 0L;
		}
		PyList_SetItem(lp, i, el_obj);
	}
	return lp;
}
/**
 * (internal)  Converts a rectangular 3D matrix of floats to a python list of
 * 	list of floats.
 *
 *	\assert(arr != 0L || len == 0)
 */
PyObject * double3DMatToList(const double ***arr, unsigned n_mats, unsigned n_rows, unsigned n_cols) {
	PyObject *lp;
	PyObject *el_obj;
	unsigned i;
	assert(arr != 0L || n_mats == 0);
	lp = PyList_New(n_mats);
	if (lp == 0L)
		return 0L;
	for (i = 0; i < n_mats; ++i) {
		el_obj = doubleMatToList(arr[i], n_rows, n_cols);
		if (el_obj == 0L) {
			Py_DECREF(lp);
			return 0L;
		}
		PyList_SetItem(lp, i, el_obj);
	}
	return lp;
}

PyObject * tupleToUnsignedArrayMaxSize(PyObject *tuple_obj, int *arr, unsigned n, unsigned *actualLen) {
	PyObject *item;
	unsigned pytuple_len = (unsigned) PyTuple_Size(tuple_obj);
	unsigned i;
	long lval;
	if (pytuple_len > n) {
	    PYTBEAGLEHON_DEBUG_PRINTF2("error in tupleToUnsignedArrayMaxSize max=%d, tuple= %d\n", n, pytuple_len);
		PyErr_SetString(PyExc_IndexError, "tuple index out of range");
		return 0L;
	}
	for (i = 0; i < n; ++i) {
		item = PyTuple_GetItem(tuple_obj, i);
		if (item == 0L) {
    		PyErr_SetString(PyExc_TypeError, "could not extract item from tuple");
			return 0L;
		}
		Py_INCREF(item);
		if (!PyInt_Check(item)) {
			Py_DECREF(item);
    		PyErr_SetString(PyExc_TypeError, "integer expected");
			return 0L;
		}
		lval = PyInt_AsLong(item);
		if (lval >= INT_MAX || lval < 0) {
			Py_DECREF(item);
    		PyErr_SetString(PyExc_TypeError, "value out of range for an unsigned integer");
	    	return 0L;
		}
		arr[i] = (int) lval;
		Py_DECREF(item);
	}
	if (actualLen)
	    *actualLen = pytuple_len;
	return none();
}

PyObject * tupleToUnsignedArray(PyObject *tuple_obj, int *arr, unsigned n) {
	PyObject *item;
	unsigned pytuple_len = (unsigned) PyTuple_Size(tuple_obj);
	unsigned i;
	long lval;
	if (pytuple_len != n) {
		PyErr_SetString(PyExc_IndexError, "tuple index out of range");
		return 0L;
	}
	for (i = 0; i < n; ++i) {
		item = PyTuple_GetItem(tuple_obj, i);
		if (item == 0L) {
    		PyErr_SetString(PyExc_TypeError, "could not extract item from tuple");
			return 0L;
		}
		Py_INCREF(item);
		if (!PyInt_Check(item)) {
			Py_DECREF(item);
    		PyErr_SetString(PyExc_TypeError, "integer expected");
			return 0L;
		}
		lval = PyInt_AsLong(item);
		if (lval >= INT_MAX || lval < 0) {
			Py_DECREF(item);
    		PyErr_SetString(PyExc_TypeError, "value out of range for an unsigned integer");
	    	return 0L;
		}
		arr[i] = (int) lval;
		Py_DECREF(item);
	}
	return none();
}

PyObject * listToUnsignedArrayMaxSize(PyObject *list_obj, int *arr, unsigned n, unsigned *actualLen) {
	PyObject *item;
	unsigned pylist_len = (unsigned) PyList_Size(list_obj);
	unsigned i;
	long lval;
	if (pylist_len > n) {
	    PYTBEAGLEHON_DEBUG_PRINTF2("error in listToUnsignedArrayMaxSize max=%d, tuple= %d\n", n, pylist_len);
		PyErr_SetString(PyExc_IndexError, "list index out of range");
		return 0L;
	}
	for (i = 0; i < n; ++i) {
		item = PyList_GetItem(list_obj, i);
		if (item == 0L) {
    		PyErr_SetString(PyExc_TypeError, "could not extract item from list");
			return 0L;
		}
		Py_INCREF(item);
		if (!PyInt_Check(item)) {
			Py_DECREF(item);
    		PyErr_SetString(PyExc_TypeError, "integer expected");
			return 0L;
		}
		lval = PyInt_AsLong(item);
		if (lval >= INT_MAX || lval < 0) {
			Py_DECREF(item);
    		PyErr_SetString(PyExc_TypeError, "value out of range for an unsigned integer");
	    	return 0L;
		}
		arr[i] = (int) lval;
		Py_DECREF(item);
	}
	if (actualLen)
	    *actualLen = pylist_len;
	return none();
}

PyObject * listToUnsignedArray(PyObject *list_obj, int *arr, unsigned n) {
	PyObject *item;
	unsigned pylist_len = (unsigned) PyList_Size(list_obj);
	unsigned i;
	long lval;
	if (pylist_len != n) {
		PyErr_SetString(PyExc_IndexError, "list index out of range");
		return 0L;
	}
	for (i = 0; i < n; ++i) {
		item = PyList_GetItem(list_obj, i);
		if (item == 0L) {
    		PyErr_SetString(PyExc_TypeError, "could not extract item from list");
			return 0L;
		}
		Py_INCREF(item);
		if (!PyInt_Check(item)) {
			Py_DECREF(item);
    		PyErr_SetString(PyExc_TypeError, "integer expected");
			return 0L;
		}
		lval = PyInt_AsLong(item);
		if (lval >= INT_MAX || lval < 0) {
			Py_DECREF(item);
    		PyErr_SetString(PyExc_TypeError, "value out of range for an unsigned integer");
	    	return 0L;
		}
		arr[i] = (int) lval;
		Py_DECREF(item);
	}
	return none();
}

PyObject * listToDoubleArrayMaxSize(PyObject *list_obj, double *arr, unsigned maxLen, unsigned *actualLen) {
	PyObject *item, *f_item;
	unsigned pylist_len = (unsigned) PyList_Size(list_obj);
	unsigned i;
	if (pylist_len > maxLen) {
		PyErr_SetString(PyExc_IndexError, "list index out of range");
		return 0L;
	}
	for (i = 0; i < maxLen; ++i) {
		item = PyList_GetItem(list_obj, i);
		if (item == 0L) {
			return 0L;
		}
		Py_INCREF(item);
		f_item = PyNumber_Float(item);
		if (f_item == 0L) {
			Py_DECREF(item);
			return 0L;
		}
		arr[i] = PyFloat_AsDouble(item);
		Py_DECREF(item);
		Py_DECREF(f_item);
	}
	if (actualLen)
	    *actualLen = pylist_len;
	return none();
}

PyObject * listToDoubleArray(PyObject *list_obj, double *arr, unsigned n) {
	PyObject *item, *f_item;
	unsigned pylist_len = (unsigned) PyList_Size(list_obj);
	unsigned i;
	if (pylist_len != n) {
		PyErr_SetString(PyExc_IndexError, "list index out of range");
		return 0L;
	}
	for (i = 0; i < n; ++i) {
		item = PyList_GetItem(list_obj, i);
		if (item == 0L) {
			return 0L;
		}
		Py_INCREF(item);
		f_item = PyNumber_Float(item);
		if (f_item == 0L) {
			Py_DECREF(item);
			return 0L;
		}
		arr[i] = PyFloat_AsDouble(item);
		Py_DECREF(item);
		Py_DECREF(f_item);
	}
	return none();
}

PyObject * listToDoubleMatrix(PyObject *list_obj, double **arr, unsigned n_rows, unsigned n_cols) {
	PyObject *item, *r_item;
	unsigned pylist_len = (unsigned) PyList_Size(list_obj);
	unsigned i;
	if (pylist_len != n_rows) {
		PyErr_SetString(PyExc_IndexError, "list index out of range");
		return 0L;
	}
	for (i = 0; i < n_rows; ++i) {
		item = PyList_GetItem(list_obj, i);
		if (item == 0L) {
			return 0L;
		}
		Py_INCREF(item);
		if (!PyList_Check(item)) {
			PyErr_SetString(PyExc_TypeError, "list of list of floats expected");
			Py_DECREF(item);
			return 0L;
		}
		r_item = listToDoubleArray(item, arr[i], n_cols);
		Py_DECREF(item);
		if (0L == r_item)
			return 0L;
	}
	return none();
}


PyObject * tupleToDoubleArray(PyObject *tuple_obj, double *arr, unsigned n, int demandExactLen) {
	PyObject *item, *f_item;
	unsigned pytuple_len = (unsigned) PyTuple_Size(tuple_obj);
	unsigned i;
	if (pytuple_len != n && ((demandExactLen != 0) || (pytuple_len < n))) {
		PyErr_SetString(PyExc_IndexError, "tuple index out of range");
		return 0L;
	}
	for (i = 0; i < n; ++i) {
		item = PyTuple_GetItem(tuple_obj, i);
		if (item == 0L) {
			return 0L;
		}
		Py_INCREF(item);
		f_item = PyNumber_Float(item);
		if (f_item == 0L) {
			Py_DECREF(item);
			return 0L;
		}
		arr[i] = PyFloat_AsDouble(item);
		Py_DECREF(item);
		Py_DECREF(f_item);
	}
	return none();
}


PyObject * tupleToDoubleMatrix(PyObject *tuple_obj, double **arr, unsigned n_rows, unsigned n_cols, int demandExactLen) {
	PyObject *item, *r_item;
	unsigned pytuple_len = (unsigned) PyTuple_Size(tuple_obj);
	unsigned i;
	if (pytuple_len != n_rows && ((demandExactLen != 0) || (pytuple_len < n_rows))) {
		PyErr_SetString(PyExc_IndexError, "tuple index out of range");
		return 0L;
	}
	for (i = 0; i < n_rows; ++i) {
		item = PyTuple_GetItem(tuple_obj, i);
		if (item == 0L) {
			return 0L;
		}
		Py_INCREF(item);
		if (!PyTuple_Check(item)) {
			PyErr_SetString(PyExc_TypeError, "tuple of tuple of floats expected");
			Py_DECREF(item);
			return 0L;
		}
		r_item = tupleToDoubleArray(item, arr[i], n_cols, demandExactLen);
		Py_DECREF(item);
		if (0L == r_item)
			return 0L;
	}
	return none();
}
