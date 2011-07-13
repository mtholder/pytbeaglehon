#include "asrv.h"
#include "phylo_util.h"
#include "py_calc_instance.h"
#include "py_asrv.h"
#include "py_util.h"
#include "calc_instance.h"
#include "py_calc_instance.h"
#include "discrete_state_model.h"
#include "internal_like_calc_env.h"

PyObject * tupleToOpCode(PyObject *tuple_obj, BeagleOperation * opPtr);

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
	PyObject * patternWeightTuple = 0L;
	double * patternWeights = 0L;
	ASRVObj ** asrvObjectArray = 0L;
	PyObject * asrvTuple = 0L;
	PyObject * fItem = 0L;
	PyObject * item = 0L;
	unsigned listSize, i;
	if (!PyArg_ParseTuple(args, "ilO!iiiiO!iiiill", &numLeaves,
											  &numPatterns,
											  &PyTuple_Type,
											  &patternWeightTuple,
											  &numStates,
											  &numStateCodeArrays,
											  &numPartialStructs,
											  &numInstRateModels,
											  &PyTuple_Type,
											  &asrvTuple,
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
    listSize = (unsigned) PyTuple_Size(patternWeightTuple);
	if (listSize > 0) {
	    if (listSize < numPatterns) {
            PyErr_SetString(PyExc_IndexError, "The length of the patternWeightTuple cannot be less than 'numPatterns'");
            return 0L;
        }
        patternWeights = (double *)malloc(numPatterns*sizeof(double));
        if (patternWeights == 0) {
    		PYTBEAGLEHON_DEBUG_PRINTF("Could not allocate patternWeights in cpytbeaglehon_init\n");
	    	PyErr_NoMemory();
		    goto errorExit;
        }
        for (i = 0; i < numPatterns; ++i) {
            item = PyTuple_GetItem(patternWeightTuple, i);
            if (item == 0L) {
                PyErr_SetString(PyExc_IndexError, "Could not extract an item from patternWeightTuple");
    		    goto errorExit;
            }
            Py_INCREF(item);
            fItem = PyNumber_Float(item);
            if (fItem == 0L) {
                Py_DECREF(item);
                PyErr_SetString(PyExc_IndexError, "Could not extract a float from patternWeightTuple");
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
    listSize = (unsigned) PyTuple_Size(asrvTuple);
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
            item = PyTuple_GetItem(asrvTuple, i);
            if (item == 0L) {
                PyErr_SetString(PyExc_IndexError, "Could not extract an item from asrvTuple");
    		    goto errorExit;
            }
            Py_INCREF(item);
            if (!PyType_IsSubtype(item->ob_type, &asrv_type)) {
                if (Py_None == (PyObject *) item ) {
                    asrvObjectArray[i] = 0L;
                    Py_DECREF(item);
                } 
                else {
                    Py_DECREF(item);
                    PyErr_SetString(PyExc_IndexError, "Could not extract a ASRVObj from asrvTuple");
    	    	    goto errorExit;
    	    	}
            }
            else {
                asrvObjectArray[i] = (ASRVObj *) item;
            }
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
	PYTBEAGLEHON_DEBUG_PRINTF("Calling createLikelihoodCalcInstance\n");
	
	
	/* \TEMP only supporting ALWAYS mode for scaling*/
	resourcePref &= (~BEAGLE_FLAG_SCALING_ALWAYS);
	resourcePref &= (~BEAGLE_FLAG_SCALING_AUTO);
    resourcePref &= (~BEAGLE_FLAG_SCALING_DYNAMIC);
	resourcePref |= BEAGLE_FLAG_SCALING_MANUAL;

	resourceReq &= (~BEAGLE_FLAG_SCALING_ALWAYS);
	resourceReq &= (~BEAGLE_FLAG_SCALING_AUTO);
    resourceReq &= (~BEAGLE_FLAG_SCALING_DYNAMIC);
	resourceReq |= BEAGLE_FLAG_SCALING_MANUAL;
	
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

#   if defined(API_TRACE_PRINTING) && API_TRACE_PRINTING
        if (handle > 0) {
            PYTBEAGLEHON_DEBUG_PRINTF("/* cAPI Call */ handle=0 is assumed in tracing code! ;\n");
        }
        if (asrvObjectArray != 0L) {
            for (i = 0; i < numInstRateModels; ++i) {
                if (asrvObjectArray[i] != 0L) {
                    PYTBEAGLEHON_DEBUG_PRINTF("/* cAPI Call */ need to create some asrv object here ;\n");
                }
            }
        }
        if (patternWeights) {
            PYTBEAGLEHON_DEBUG_PRINTF("/* cAPI Call */ need to create some pattern weights here ;\n");
        }
    	PYTBEAGLEHON_DEBUG_PRINTF3("/* cAPI Call */ int numLeaves=%d; long numPatterns=%ld ; const double * patternWeights=0L; int numStates=%d;",  numLeaves, numPatterns, numStates);
	    PYTBEAGLEHON_DEBUG_PRINTF3("int numStateCodeArrays=%d; int numPartialStructs=%d; int numInstRateModels=%d ; const ASRVObj ** asrvObjectArray=0L;",  numStateCodeArrays, numPartialStructs, numInstRateModels);
    	PYTBEAGLEHON_DEBUG_PRINTF4("int numProbMats=%d; int numEigenStorage=%d; int numRescalingsMultipliers=%d; int resourceIndex=%d; ",  numProbMats, numEigenStorage, numRescalingsMultipliers, resourceIndex);
	    PYTBEAGLEHON_DEBUG_PRINTF2("long resourcePref=%ld; long resourceReq=%ld;\n", resourcePref, resourceReq);
    	PYTBEAGLEHON_DEBUG_PRINTF("/* cAPI Call */ long handle = createLikelihoodCalcInstance(numLeaves, numPatterns, patternWeights, numStates, numStateCodeArrays, numPartialStructs, numInstRateModels, asrvObjectArray, numProbMats, numEigenStorage, numRescalingsMultipliers, resourceIndex, resourcePref, resourceReq);\n");
    	PYTBEAGLEHON_DEBUG_PRINTF1("/* cAPI Call */ const struct LikeCalculatorInstance * LCI = getLikeCalculatorInstance(%ld); if (LCI == 0L) {fprintf(stderr, \"getLikeCalculatorInstance failed\"); return 1;}\n", handle);
#   endif

	if (handle < 0) {
		PYTBEAGLEHON_DEBUG_PRINTF("Could not allocate likelihood_calculator_instance\n");
		PyErr_NoMemory();
		goto errorExit;
	}
	PYTBEAGLEHON_DEBUG_PRINTF1("Returning handle for instance %ld\n", handle);
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


PyObject* cPytBeagleHonFree(PyObject *self, PyObject *args) {
	long handle;
	if (!PyArg_ParseTuple(args, "l", &handle)) {
		return 0L;
	}
#   if defined(API_TRACE_PRINTING) && API_TRACE_PRINTING
        PYTBEAGLEHON_DEBUG_PRINTF("/* cAPI Call */ freeLikeCalculatorInstance(handle);\n");
#   endif
	if (freeLikeCalculatorInstance(handle) != 0) {
		PyErr_SetString(PyExc_ValueError, "Error freeing calculation instance");
		return 0L;
	}
	return none();
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

PyObject* pyGetModelList(PyObject *self, PyObject *args) {
	long handle;
	unsigned numModels, i;
	PyObject * toReturn;
    const DSCTModelObj ** modArray;
    if (!PyArg_ParseTuple(args, "l", &handle))
		return 0L;
    modArray = getModelList(handle, &numModels);
    if (modArray == 0L) {
        PyErr_SetString(PyExc_IndexError, "Invalid likelihood instance index");
        return 0L;
    }
	toReturn = PyTuple_New(numModels);
	for (i = 0; i < numModels; ++i) {
	    PyTuple_SetItem(toReturn, i, (PyObject *)(modArray[i]));
	}
    return toReturn;
}

PyObject* pySetStateCodeArray(PyObject *self, PyObject *args) {
#   if defined(DEBUG_PRINTING) && DEBUG_PRINTING
        int i;        
#   endif
	long handle;
	int stateCodeArrayIndex;
	unsigned scArraySize;
	unsigned numCopied;
    PyObject * scTuple = 0L;
    if (!PyArg_ParseTuple(args, "liO!", &handle, &stateCodeArrayIndex, &PyTuple_Type, &scTuple))
        return 0L;
    scArraySize = (unsigned) PyTuple_Size(scTuple);
    const struct LikeCalculatorInstance * LCI = getLikeCalculatorInstance(handle);
    if (LCI == 0L) {
        PyErr_SetString(PyExc_IndexError, "Invalid likelihood instance index");
        return 0L;
    }
    if (scArraySize != LCI->numPatterns) {
        PyErr_SetString(PyExc_IndexError, "Length of state code array passed in exceeds capacity of the likelihood calculator instances");
        return 0L;
    }
    if (!tupleToUnsignedArrayMaxSize(scTuple, LCI->stateCodeArrayScratch, LCI->numPatterns, &numCopied))
        return 0L;

#   if defined(API_TRACE_PRINTING) && API_TRACE_PRINTING
        PYTBEAGLEHON_DEBUG_PRINTF2("/* cAPI Call */ handle = %ld; int scarr%d[] = {", handle, stateCodeArrayIndex);
        PYTBEAGLEHON_DEBUG_PRINTF1("%d", LCI->stateCodeArrayScratch[0]);
        for (i = 1; i < LCI->numPatterns; ++i) {
            PYTBEAGLEHON_DEBUG_PRINTF1(", %d", LCI->stateCodeArrayScratch[i]);
        }
        PYTBEAGLEHON_DEBUG_PRINTF("};\n");
        PYTBEAGLEHON_DEBUG_PRINTF2("/* cAPI Call */ rc = setStateCodeArray(handle, %d, scarr%d); if (rc != 0L) {fprintf(stderr, \"setStateCodeArray failed\"); return rc;}\n", stateCodeArrayIndex, stateCodeArrayIndex);
#   endif
        
    if (setStateCodeArray(handle, stateCodeArrayIndex, LCI->stateCodeArrayScratch) != 0) {
        PyErr_SetString(PyExc_IndexError, "Error calling setStateCodeArray");
        return 0L;
    }
    return none();


}

PyObject* pyCalcPartials(PyObject *self, PyObject *args) {
#   if defined(DEBUG_PRINTING) && DEBUG_PRINTING
        BeagleOperation * bo;
#   endif
	long handle;
	unsigned i;
	unsigned opArraySize, waitTupleSize;
	PyObject *item, *r_item;
    PyObject * opList = 0L;
    PyObject * waitTuple = 0L;
    if (!PyArg_ParseTuple(args, "lO!O!", &handle, &PyList_Type, &opList, &PyTuple_Type, &waitTuple))
        return 0L;
    opArraySize = (unsigned) PyList_Size(opList);
    waitTupleSize = (unsigned) PyTuple_Size(waitTuple);
    const struct LikeCalculatorInstance * LCI = getLikeCalculatorInstance(handle);
    if (LCI == 0L) {
        PyErr_SetString(PyExc_IndexError, "Invalid likelihood instance index");
        return 0L;
    }
    if (opArraySize > LCI->numPartialStructs) {
        PyErr_SetString(PyExc_IndexError, "Length number of update partial operations in one call cannot exceed the number of partial structures requested for the instance.");
        return 0L;
    }
    if (waitTupleSize > LCI->numPartialStructs) {
        PYTBEAGLEHON_DEBUG_PRINTF2("waitTupleSize=%d LCI->numPartialStructs=%d\n", waitTupleSize, LCI->numPartialStructs);
        PyErr_SetString(PyExc_IndexError, "The number of partials to wait for cannot exceed the number of partial structures requested for the instance.");
        return 0L;
    }

	for (i = 0; i < opArraySize; ++i) {
		item = PyList_GetItem(opList, i);
		if (item == 0L) {
			return 0L;
		}
		Py_INCREF(item);
		if (!PyTuple_Check(item)) {
			PyErr_SetString(PyExc_TypeError, "list of operation list expected");
			Py_DECREF(item);
			return 0L;
		}
		r_item = tupleToOpCode(item, LCI->opScratch + i);
		Py_DECREF(item);
		if (0L == r_item)
			return 0L;
	}
	tupleToUnsignedArrayMaxSize(waitTuple, LCI->waitPartialIndexScratch, LCI->numPartialStructs, &waitTupleSize);

#   if defined(API_TRACE_PRINTING) && API_TRACE_PRINTING
        PYTBEAGLEHON_DEBUG_PRINTF1("/* cAPI Call */ handle = %ld; ", handle);
        for (i = 0; i < opArraySize; ++i) {
            bo = LCI->opScratch + i;
            PYTBEAGLEHON_DEBUG_PRINTF4("LCI->opScratch[%d] = partialOperation(%d, %d, %d, ", i,  bo->destinationPartials, bo->destinationScaleWrite, bo->destinationScaleRead);
            PYTBEAGLEHON_DEBUG_PRINTF4("%d, %d, %d, %d); ", bo->child1Partials, bo->child1TransitionMatrix, bo->child2Partials, bo->child2TransitionMatrix);
        }
        PYTBEAGLEHON_DEBUG_PRINTF("\n/* cAPI Call */ ");
        if (waitTupleSize < 1) {
            PYTBEAGLEHON_DEBUG_PRINTF1("/* cAPI Call */ rc = calcPartials(handle, LCI->opScratch, %d, 0, -1);\n", opArraySize);
        }
        else {
            for (i = 0; i < waitTupleSize; ++i) {
                PYTBEAGLEHON_DEBUG_PRINTF2("LCI->waitPartialIndexScratch[%d] = %d; ", i, LCI->waitPartialIndexScratch[i]);
            }
            PYTBEAGLEHON_DEBUG_PRINTF2("\n/* cAPI Call */ rc = calcPartials(handle, LCI->opScratch, %d, LCI->waitPartialIndexScratch, %d);\n", opArraySize, waitTupleSize);
        }
        PYTBEAGLEHON_DEBUG_PRINTF("/* cAPI Call */ if (rc != 0L) {fprintf(stderr, \"calcPartials failed\"); return rc;}\n");
#   endif


	if (calcPartials(handle, LCI->opScratch, opArraySize, LCI->waitPartialIndexScratch, (int) waitTupleSize) != 0) {
        PyErr_SetString(PyExc_ValueError, "Error calling calcPartials");
	    return 0L;
	}
	return none();
}


PyObject * tupleToOpCode(PyObject *tuple_obj, BeagleOperation * opPtr) {
    assert(opPtr);
	unsigned pytuple_len = (unsigned) PyTuple_Size(tuple_obj);
	long lval;
	int tmp;
	if (pytuple_len != 7) {
		PyErr_SetString(PyExc_IndexError, "Operation tuple should have 7 integers");
		return 0L;
	}
    /* dest partial */
    if (extractNonNegativeIntFromTuple(tuple_obj, 0, &tmp) == 0)
        return 0L;
    opPtr->destinationPartials = tmp;
    if (extractLongFromTuple(tuple_obj, 1, &lval) == 0)
        return 0L;
	opPtr->destinationScaleWrite = (int) lval;
    if (extractLongFromTuple(tuple_obj, 2, &lval) == 0)
        return 0L;
	opPtr->destinationScaleRead = lval;
    if (extractNonNegativeIntFromTuple(tuple_obj, 3, &tmp) == 0)
        return 0L;
	opPtr->child1Partials = tmp;
    if (extractNonNegativeIntFromTuple(tuple_obj, 4, &tmp) == 0)
        return 0L;
	opPtr->child1TransitionMatrix = tmp;
    if (extractNonNegativeIntFromTuple(tuple_obj, 5, &tmp) == 0)
        return 0L;
	opPtr->child2Partials = tmp;
    if (extractNonNegativeIntFromTuple(tuple_obj, 6, &tmp) == 0)
        return 0L;
	opPtr->child2TransitionMatrix = tmp;
	return none();
}


PyObject* pySetStateFreq(PyObject *self, PyObject *args) {
	long handle;
	unsigned eigenInd;
	unsigned weightTupleSize;
    PyObject * weightTuple = 0L;
    if (!PyArg_ParseTuple(args, "liO!", &handle, &eigenInd, &PyTuple_Type, &weightTuple))
        return 0L;
    weightTupleSize = (unsigned) PyTuple_Size(weightTuple);
    const struct LikeCalculatorInstance * LCI = getLikeCalculatorInstance(handle);
    if (LCI == 0L) {
        PyErr_SetString(PyExc_IndexError, "Invalid likelihood instance index");
        return 0L;
    }
    if (weightTupleSize != LCI->numStates) {
        PyErr_SetString(PyExc_IndexError, "The length of the state frequency array does not equal the number of states");
        return 0L;
    }
    if (eigenInd <0 || eigenInd >= LCI->numEigenStorage) {
        PyErr_SetString(PyExc_IndexError, "The index of the buffer for the state frequencies is out of range.");
        return 0L;
    }
    if (tupleToDoubleArray(weightTuple, LCI->categWeightScratch, weightTupleSize, 1) == 0L)
        return 0L;
#   if defined(API_TRACE_PRINTING) && API_TRACE_PRINTING
        int i;
        PYTBEAGLEHON_DEBUG_PRINTF2("/* cAPI Call */ handle = %ld; eigenIndex = %d; ", handle, eigenInd);
        for (i = 0; i < weightTupleSize; ++i) {
            PYTBEAGLEHON_DEBUG_PRINTF2("LCI->categWeightScratch[%d] = %.10lf; ", i,  LCI->categWeightScratch[i]);
        }
        PYTBEAGLEHON_DEBUG_PRINTF("\n/* cAPI Call */ rc = setStateFreq(handle, eigenIndex, LCI->categWeightScratch);\n");
        PYTBEAGLEHON_DEBUG_PRINTF("/* cAPI Call */ if (rc != 0L) {fprintf(stderr, \"setStateFreq failed\"); return rc;}\n");
#   endif
	if (setStateFreq(handle, eigenInd, LCI->categWeightScratch) != 0) {
        PyErr_SetString(PyExc_ValueError, "Error calling setStateFreq");
	    return 0L;
	}
	return none();
}



PyObject* pyCalcRootLnLikelihood(PyObject *self, PyObject *args) {
	long handle;
	unsigned indArraySize, i;
	double lnL;
    PyObject * rootPartialIndTuple = 0L;
    PyObject * categWeightIndTuple = 0L;
    PyObject * stateFreqIndTuple = 0L;
    PyObject * rescalerIndTuple = 0L;
    if (!PyArg_ParseTuple(args, "lO!O!O!O!", &handle, &PyTuple_Type, &rootPartialIndTuple, 
                                                  &PyTuple_Type, &categWeightIndTuple, 
                                                  &PyTuple_Type, &stateFreqIndTuple, 
                                                  &PyTuple_Type, &rescalerIndTuple))
        return 0L;
    indArraySize = (unsigned) PyTuple_Size(rootPartialIndTuple);
    i = (unsigned) PyTuple_Size(categWeightIndTuple);
    if (indArraySize != i){
        PyErr_SetString(PyExc_IndexError, "Expecting the number of root partials to equal the number of category weight buffer indices");
        return 0L;
    }
    i = (unsigned) PyTuple_Size(stateFreqIndTuple);
    if (indArraySize != i){
        PyErr_SetString(PyExc_IndexError, "Expecting the number of root partials to equal the number of state frequency buffer indices");
        return 0L;
    }
    i = (unsigned) PyTuple_Size(rescalerIndTuple);
    if (indArraySize != i){
        PyErr_SetString(PyExc_IndexError, "Expecting the number of root partials to equal the number of rescaler buffer indices");
        return 0L;
    }
    const struct LikeCalculatorInstance * LCI = getLikeCalculatorInstance(handle);
    if (LCI == 0L) {
        PyErr_SetString(PyExc_IndexError, "Invalid likelihood instance index");
        return 0L;
    }
    if (indArraySize > LCI->numEigenStorage) {
        PyErr_SetString(PyExc_IndexError, "Oddly enough, the number of rate category weights has to be less than the number of eigen solution storage slots.");
        return 0L;
    }
    if (tupleToUnsignedArray(rootPartialIndTuple, LCI->rootPartialIndexScratch, indArraySize) == 0L)
        return 0L;
    if (tupleToUnsignedArray(categWeightIndTuple, LCI->categWeightIndexScratch, indArraySize) == 0L)
        return 0L;
    if (tupleToUnsignedArray(stateFreqIndTuple, LCI->stateFreqIndexScratch, indArraySize) == 0L)
        return 0L;
    if (tupleToIntArray(rescalerIndTuple, LCI->rootRescalerIndexScratch, indArraySize) == 0L)
        return 0L;
#   if defined(API_TRACE_PRINTING) && API_TRACE_PRINTING
        PYTBEAGLEHON_DEBUG_PRINTF("/* cAPI Call */ ");
        for (i = 0; i < indArraySize; ++i) {
            PYTBEAGLEHON_DEBUG_PRINTF2("LCI->rootPartialIndexScratch[%d] = %d; ", i,  LCI->rootPartialIndexScratch[i]);
        }
        PYTBEAGLEHON_DEBUG_PRINTF("/* cAPI Call */ ");
        for (i = 0; i < indArraySize; ++i) {
            PYTBEAGLEHON_DEBUG_PRINTF2("LCI->categWeightIndexScratch[%d] = %d; ", i,  LCI->categWeightIndexScratch[i]);
        }
        PYTBEAGLEHON_DEBUG_PRINTF("/* cAPI Call */ ");
        for (i = 0; i < indArraySize; ++i) {
            PYTBEAGLEHON_DEBUG_PRINTF2("LCI->stateFreqIndexScratch[%d] = %d; ", i,  LCI->stateFreqIndexScratch[i]);
        }
        PYTBEAGLEHON_DEBUG_PRINTF("/* cAPI Call */ ");
        for (i = 0; i < indArraySize; ++i) {
            PYTBEAGLEHON_DEBUG_PRINTF2("LCI->rootRescalerIndexScratch[%d] = %d; ", i,  LCI->rootRescalerIndexScratch[i]);
        }
        PYTBEAGLEHON_DEBUG_PRINTF1("\n/* cAPI Call */ rc = calcRootLnL(handle, LCI->rootPartialIndexScratch, LCI->categWeightIndexScratch, LCI->stateFreqIndexScratch, LCI->rootRescalerIndexScratch, %d, &lnL);\n", indArraySize);
        PYTBEAGLEHON_DEBUG_PRINTF("/* cAPI Call */ if (rc != 0L) {fprintf(stderr, \"calcRootLnL failed\"); return rc;}\n");
#   endif
	if (calcRootLnL(handle, LCI->rootPartialIndexScratch, LCI->categWeightIndexScratch, LCI->stateFreqIndexScratch, LCI->rootRescalerIndexScratch, (int) indArraySize, &lnL) != 0) {
        PyErr_SetString(PyExc_ValueError, "Error calling calcRootLnL");
	    return 0L;
	}
	PYTBEAGLEHON_DEBUG_PRINTF1("calcRootLnL => %.10lf\n", lnL);
	return PyFloat_FromDouble(lnL);
}


PyObject* pySetSingletonCatWts(PyObject *self, PyObject *args) {
	long handle;
	unsigned indArraySize, weightTupleSize;
    PyObject * indList = 0L;
    PyObject * weightTuple = 0L;
    if (!PyArg_ParseTuple(args, "lO!O!", &handle, &PyTuple_Type, &indList, &PyTuple_Type, &weightTuple))
        return 0L;
    indArraySize = (unsigned) PyTuple_Size(indList);
    weightTupleSize = (unsigned) PyTuple_Size(weightTuple);
    if (weightTupleSize != indArraySize){
        PyErr_SetString(PyExc_IndexError, "Expecting the index list and weight list to have the same length");
        return 0L;
    }
    const struct LikeCalculatorInstance * LCI = getLikeCalculatorInstance(handle);
    if (LCI == 0L) {
        PyErr_SetString(PyExc_IndexError, "Invalid likelihood instance index");
        return 0L;
    }
    if (indArraySize > LCI->numEigenStorage) {
        PyErr_SetString(PyExc_IndexError, "Oddly enough, the number of rate category weights has to be less than the number of eigen solution storage slots.");
        return 0L;
    }
    if (tupleToUnsignedArray(indList, LCI->categWeightIndexScratch, indArraySize) == 0L)
        return 0L;
    if (tupleToDoubleArray(weightTuple, LCI->categWeightScratch, indArraySize, 1) == 0L)
        return 0L;
#   if defined(API_TRACE_PRINTING) && API_TRACE_PRINTING
        int i;
        PYTBEAGLEHON_DEBUG_PRINTF1("/* cAPI Call */ handle = %ld; ", handle);
        for (i = 0; i < indArraySize; ++i) {
            PYTBEAGLEHON_DEBUG_PRINTF2("LCI->categWeightIndexScratch[%d] = %d; ", i,  LCI->categWeightIndexScratch[i]);
        }
        for (i = 0; i < indArraySize; ++i) {
            PYTBEAGLEHON_DEBUG_PRINTF2("LCI->categWeightScratch[%d] = %.10lf; ", i,  LCI->categWeightScratch[i]);
        }
        PYTBEAGLEHON_DEBUG_PRINTF1("\n/* cAPI Call */ rc = setSingletonCategoryWeights(handle, LCI->categWeightIndexScratch, LCI->categWeightScratch, %d);\n", (int) indArraySize);
        PYTBEAGLEHON_DEBUG_PRINTF("/* cAPI Call */ if (rc != 0L) {fprintf(stderr, \"setSingletonCategoryWeights failed\"); return rc;}\n");
#   endif
	if (setSingletonCategoryWeights(handle, LCI->categWeightIndexScratch, LCI->categWeightScratch, (int) indArraySize) != 0) {
        PyErr_SetString(PyExc_ValueError, "Error calling setSingletonCategoryWeights");
	    return 0L;
	}
	return none();
}

