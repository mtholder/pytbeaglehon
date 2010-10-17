/**  Python-interop hooks for among-site rate variation stuff.
 */
#if ! defined(PY_DISCRETE_STATE_MODEL_H)
#define PY_DISCRETE_STATE_MODEL_H

#include "pytbeaglehon_defs.h"

#ifdef __cplusplus
extern "C" 
{
#endif

extern PyTypeObject dsct_model_type;

PyObject* cdsctm_set_q_mat(PyObject *self, PyObject *args);
PyObject* cdsctm_calc_eigens(PyObject *self, PyObject *args);
PyObject* cdsctm_calc_pr_mats(PyObject *self, PyObject *args);


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

