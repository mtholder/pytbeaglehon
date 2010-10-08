/**
  adaptor code for filling in ASRVObj

    This file does not actually use python-c extension (calls like PyObject_Del are
        defined to free in pytbeaglehon_defs.h if BUILDING_FOR_PYTHON is not defined).
*/
#include <stdlib.h>
#include "asrv.h"
#include "phylo_util.h"

void internal_asrv_set_shape(ASRVObj *asrh, double alpha) {
	PYTBEAGLEHON_DEBUG_PRINTF("In internal_asrv_set_shape\n");
	PYTBEAGLE_ASSERT(asrh);
	PYTBEAGLE_ASSERT(alpha > 0.0);
	double beta;
	asrh->param = alpha;
	if (asrh->n < 1)
	    return;
	beta = 1.0/alpha;
	DiscreteGamma(asrh->freq, asrh->rate, alpha, beta, asrh->n, asrh->style);
}



void asrv_obj_dtor(ASRVObj * asrh) {
	PYTBEAGLEHON_DEBUG_PRINTF("In asrh dtor\n");
	if (asrh == 0L)
		return;
	if (asrh->rate != 0L)
		free(asrh->rate);
	if (asrh->freq != 0L)
		free(asrh->freq);
	PyObject_Del(asrh);
}




/* fills in the fields of asrh and returns it (or returns null if there was
    an allocation error filling the fields).
*/
ASRVObj * private_asrv_obj_init(ASRVObj * asrh, unsigned dim, int style, double param) {
    PYTBEAGLE_ASSERT(asrh);
    const unsigned arr_len = dim*sizeof(double);
    asrh->n = dim;
    asrh->rate = 0L;
    asrh->freq = 0L;
    asrh->style = style;
    asrh->param = param;
    if (dim > 0) {
        asrh->rate = (double*)malloc(arr_len*sizeof(double));
        if (asrh->rate == 0L) {
            PyErr_NoMemory();
            return 0L;
        }
        asrh->freq = (double*)malloc(arr_len*sizeof(double));
        if (asrh->freq == 0L) {
            return 0L;
        }
    internal_asrv_set_shape(asrh, param);
    }
    return asrh;
}





#if defined(BUILDING_FOR_PYTHON) && BUILDING_FOR_PYTHON
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

#endif  /*defined(BUILDING_FOR_PYTHON) && BUILDING_FOR_PYTHON */

ASRVObj* asrv_obj_new(unsigned dim, int style, double param)  {
    PYTBEAGLE_ASSERT(dim > 1);
    PYTBEAGLEHON_DEBUG_PRINTF("In asrv_obj_new\n");
    ASRVObj * asrh = PyObject_New(ASRVObj, &asrv_type);
    if (asrh) {
        if (private_asrv_obj_init(asrh, dim, style, param) == 0L)
            goto errorExit;
    }
    else {
        PYTBEAGLEHON_DEBUG_PRINTF("PyObject_New returned 0L\n");
    }
    return asrh;
    errorExit:
        PYTBEAGLEHON_DEBUG_PRINTF("In asrv_obj_new errorExit\n");
        asrv_obj_dtor(asrh);
        return 0L;
}





/**
  pytbeaglehon phylogenetic likelihood caluclations using beaglelib.

  Copyright 2010 Mark T. Holder (except where noted). All rights reserved.

  See "LICENSE.txt" for terms and conditions of usage.
*/
