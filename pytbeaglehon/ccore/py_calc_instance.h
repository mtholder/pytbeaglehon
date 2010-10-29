/**  Python-interop hooks for configuring a "calculation instance" that holds
    the data structures needed for likelihoo calculations.
 */
#if ! defined(PY_CALC_INSTANCE_H)
#define PY_CALC_INSTANCE_H

#include "pytbeaglehon_defs.h"

#ifdef __cplusplus
extern "C" 
{
#endif

PyObject* cPytBeagleHonFree(PyObject *self, PyObject *args);
PyObject* cPytBeagleHonInit(PyObject *self, PyObject *args);
PyObject* pyGetNumComputationalResources(PyObject *self, PyObject *args);
PyObject* pyGetResourceInfo(PyObject *self, PyObject *args);
PyObject* pyGetModelList(PyObject *self, PyObject *args);
PyObject* pySetStateCodeArray(PyObject *self, PyObject *args);
PyObject* pyCalcPartials(PyObject *self, PyObject *args);
PyObject* pySetSingletonCatWts(PyObject *self, PyObject *args);
PyObject* pySetStateFreq(PyObject *self, PyObject *args);
PyObject* pyCalcRootLnLikelihood(PyObject *self, PyObject *args);

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

