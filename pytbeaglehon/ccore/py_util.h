/**  generic helper functions for writing python-c extensions
 */
#if ! defined(PY_UTIL_H)
#define PY_UTIL_H

#include "pytbeaglehon_defs.h"

#ifdef __cplusplus
extern "C" 
{
#endif



/*******************************************************************************
 * Function definitions
 */

		/***********************************************************************
		 * Utility functions for Python binding
		 */

PyObject* none(void);

/**
 * (internal) converts an python index on [-dim,dim) to the positive version
 * of the index.
 * Returns a non-negative index for success or -1.
 * if -1 is returned then the function will have generated a
 * IndexError("list index out of range") already.
 *
 */
int pyIndToCInd(int py_ind, unsigned dim);
/**
 * (internal) Converts a double array to a python list of floats
 *
 *	\assert(arr != 0L || len == 0)
 */
PyObject * doubleArrayToList(const double *arr, unsigned len);
PyObject * doubleMatToList(const double **arr, unsigned n_rows, unsigned n_cols);
PyObject * double3DMatToList(const double ***arr, unsigned n_mats, unsigned n_rows, unsigned n_cols);
PyObject * tupleToUnsignedArrayMaxSize(PyObject *list_obj, int *arr, unsigned n, unsigned *actualLen);
PyObject * tupleToUnsignedArray(PyObject *list_obj, int *arr, unsigned n);
PyObject * listToUnsignedArrayMaxSize(PyObject *list_obj, int *arr, unsigned n, unsigned *actualLen);
PyObject * listToUnsignedArray(PyObject *list_obj, int *arr, unsigned n);
PyObject * listToDoubleArrayMaxSize(PyObject *list_obj, double *arr, unsigned maxLen, unsigned *actualLen);
PyObject * listToDoubleArray(PyObject *list_obj, double *arr, unsigned n);
PyObject * listToDoubleMatrix(PyObject *list_obj, double **arr, unsigned n_rows, unsigned n_cols);
PyObject * tupleToDoubleArray(PyObject *tuple_obj, double *arr, unsigned n, int demandExactLen);
PyObject * tupleToDoubleMatrix(PyObject *tuple_obj, double **arr, unsigned n_rows, unsigned n_cols, int demandExactLen);
int extractLongFromTuple(PyObject * tuple_obj, unsigned i, long *value);
int extractNonNegativeIntFromTuple(PyObject * tuple_obj, unsigned i, int *value);

#ifdef __cplusplus
}
/* end of extern c bit */
#endif






#endif

/**
  pytbeaglehon phylogenetic likelihood caluclations using beaglelib.

  Copyright 2010 Mark T. Holder (except where noted). All rights reserved.

  See "LICENSE.txt" for terms and conditions of usage.
  
  
  
  Code for writing extensions was adapted from Alex Martelli's example code 
    in the Python Cookbook 17.1

*/

