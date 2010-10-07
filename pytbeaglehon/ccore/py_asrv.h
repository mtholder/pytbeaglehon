/**  Python-interop hooks for among-site rate variation stuff.
 */
#if ! defined(PY_ASRV_H)
#define PY_ASRV_H

#include "pytbeaglehon_defs.h"

#ifdef __cplusplus
extern "C" 
{
#endif

extern PyTypeObject asrv_type;

PyObject* casrvo_ctor(PyObject *self, PyObject *args);
PyObject* casrvo_get_n_cat(PyObject *self, PyObject *args);
PyObject* casrvo_get_rates(PyObject *self, PyObject *args);
PyObject* casrvo_get_shape(PyObject *self, PyObject *args);
PyObject* casrvo_set_shape(PyObject *self, PyObject *args);



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

