#include "asrv.h"
#include "phylo_util.h"
#include "py_calc_instance.h"
#include "py_asrv.h"
#include "py_util.h"
#include "calc_instance.h"
#include "py_calc_instance.h"



PyObject* cPytBeagleHonFree(PyObject *self, PyObject *args) {
	long handle;
	if (!PyArg_ParseTuple(args, "l", &handle)) {
		return 0L;
	}
	if (freeLikeCalculatorInstance(handle) != 0) {
		PyErr_SetString(PyExc_ValueError, "Error freeing calculation instance");
		return 0L;
	}
	return none();
}
PyObject* cPytBeagleHonInit(PyObject *self, PyObject *args) {
	int numLeaves;
	long numPatterns;
	int numStates;
	int numStateCodeArrays;
	int numPartialStructs;
	int numInstRateModels;
	int numProbMats;
	int numEigenStorage;
	int numRescalingsMultipliers;
	int resourceIndex;
	long resourcePref, resourceReq;
	long handle;
	PyObject * patternWeightList = 0L;
	double * patternWeights = 0L;
	ASRVObj ** asrvObjectArray = 0L;
	PyObject * asrvList = 0L;
	PyObject * fItem = 0L;
	PyObject * item = 0L;
	unsigned listSize, i;
	if (!PyArg_ParseTuple(args, "ilO!iiiiO!iiiill", &numLeaves,
											  &numPatterns,
											  &PyList_Type,
											  &patternWeightList,
											  &numStates,
											  &numStateCodeArrays,
											  &numPartialStructs,
											  &numInstRateModels,
											  &PyList_Type,
											  &asrvList,
											  &numProbMats,
											  &numEigenStorage,
											  &numRescalingsMultipliers,
											  &resourceIndex,
											  &resourcePref,
											  &resourceReq)) {
		return 0L;
	}
	if (numLeaves < 0) {
		PyErr_SetString(PyExc_ValueError, "The number of leaves cannot be negative");
		return 0;
	}
	if (numPatterns < 0) {
		PyErr_SetString(PyExc_ValueError, "The number of patterns cannot be negative");
		return 0;
	}
    listSize = (unsigned) PyList_Size(patternWeightList);
	if (listSize > 0) {
	    if (listSize < numPatterns) {
            PyErr_SetString(PyExc_IndexError, "The length of the patternWeightList cannot be less than 'numPatterns'");
            return 0L;
        }
        patternWeights = (double *)malloc(numPatterns*sizeof(double));
        if (patternWeights == 0) {
    		PYTBEAGLEHON_DEBUG_PRINTF("Could not allocate patternWeights in cpytbeaglehon_init\n");
	    	PyErr_NoMemory();
		    goto errorExit;
        }
        for (i = 0; i < numPatterns; ++i) {
            item = PyList_GetItem(patternWeightList, i);
            if (item == 0L) {
                PyErr_SetString(PyExc_IndexError, "Could not extract an item from patternWeightList");
    		    goto errorExit;
            }
            Py_INCREF(item);
            fItem = PyNumber_Float(item);
            if (fItem == 0L) {
                Py_DECREF(item);
                PyErr_SetString(PyExc_IndexError, "Could not extract a float from patternWeightList");
    		    goto errorExit;
            }
            patternWeights[i] = PyFloat_AsDouble(item);
            Py_DECREF(item);
            Py_DECREF(fItem);
        }
	}


	if (numStates < 2) {
		PyErr_SetString(PyExc_ValueError, "The number of states cannot be less than 1");
        goto errorExit;
	}
	if (numStateCodeArrays < 0) {
		PyErr_SetString(PyExc_ValueError, "The number of state code arrays cannot be negative");
        goto errorExit;
	}
	if (numPartialStructs < 0) {
		PyErr_SetString(PyExc_ValueError, "The number of partial likelihood arrays cannot be negative");
        goto errorExit;
	}
	if (numInstRateModels < 0) {
		PyErr_SetString(PyExc_ValueError, "The number of model instantaneous relative rate matrices cannot be negative");
        goto errorExit;
	}
    listSize = (unsigned) PyList_Size(asrvList);
	if (listSize > 0) {
	    if (listSize < numInstRateModels) {
            PyErr_SetString(PyExc_IndexError, "The length of the asrv objects cannot be less than 'asrvObjectArray'");
            return 0L;
        }
        asrvObjectArray = (ASRVObj **)mallocZeroedPointerArray(numInstRateModels);
        if (asrvObjectArray == 0) {
    		PYTBEAGLEHON_DEBUG_PRINTF("Could not allocate asrvObjectArray in cpytbeaglehon_init\n");
	    	PyErr_NoMemory();
		    goto errorExit;
        }
        for (i = 0; i < numInstRateModels; ++i) {
            item = PyList_GetItem(asrvList, i);
            if (item == 0L) {
                PyErr_SetString(PyExc_IndexError, "Could not extract an item from asrvList");
    		    goto errorExit;
            }
            Py_INCREF(item);
            if (!PyType_IsSubtype(item->ob_type, &asrv_type)) {
                Py_DECREF(item);
                PyErr_SetString(PyExc_IndexError, "Could not extract a ASRVObj from asrvList");
    		    goto errorExit;
            }
            asrvObjectArray[i] = (ASRVObj *) item;
        }
	}
	if (numProbMats < 0) {
		PyErr_SetString(PyExc_ValueError, "The number of transition probability matrices cannot be negative");
        goto errorExit;
	}
	if (numEigenStorage < 0) {
		PyErr_SetString(PyExc_ValueError, "The number of eigen solutions cannot be negative");
        goto errorExit;
	}
	if (numRescalingsMultipliers < 0) {
		PyErr_SetString(PyExc_ValueError, "The number of rescaling buffers cannot be negative");
        goto errorExit;
	}
	handle = createLikelihoodCalcInstance(
	        numLeaves,
            numPatterns,
            patternWeights,
            numStates,
            numStateCodeArrays,
            numPartialStructs,
            numInstRateModels,
            (const ASRVObj **) asrvObjectArray,
            numProbMats,
            numEigenStorage,
            numRescalingsMultipliers,
            resourceIndex,
            resourcePref,
            resourceReq);
	if (handle < 0) {
		PYTBEAGLEHON_DEBUG_PRINTF("Could not allocate likelihood_calculator_instance\n");
		PyErr_NoMemory();
		goto errorExit;
	}
	return PyInt_FromLong(handle);
	errorExit:
	    if (patternWeights != 0)
	        free(patternWeights);
	    if (asrvObjectArray != 0) {
	        for (i = 0; i < numInstRateModels; ++i) {
	            if (asrvObjectArray[i]) {
	                Py_DECREF((PyObject *)asrvObjectArray[i]);
	            }
	        }
	        free(asrvObjectArray);
	    }
	    return 0L;
}


PyObject* pyGetNumComputationalResources(PyObject *self, PyObject *args) {
	unsigned numResources = getNumComputationalResources();
    long ln = (long) numResources;
    return PyInt_FromLong(ln);
}

PyObject* pyGetResourceInfo(PyObject *self, PyObject *args) {
	int resourceIndex, rc;
    char name[81];
    char description[81];
    long optsFlags, reqFlags;
	PyObject * toReturn;
	if (!PyArg_ParseTuple(args, "i", &resourceIndex))
		return 0L;
    rc = getComputationalResourceDetails(resourceIndex, name, description, &optsFlags, &reqFlags);
    if (rc < 0) {
        PyErr_SetString(PyExc_IndexError, "call to getComputationalResourceDetails failed");
        return 0L;
    }
	toReturn = PyTuple_New(4);

    PyTuple_SetItem(toReturn, 0, PyString_FromString(name));
    PyTuple_SetItem(toReturn, 1, PyString_FromString(description));
    PyTuple_SetItem(toReturn, 2, PyInt_FromLong(optsFlags));
    PyTuple_SetItem(toReturn, 3, PyInt_FromLong(reqFlags));
	return toReturn;
}

