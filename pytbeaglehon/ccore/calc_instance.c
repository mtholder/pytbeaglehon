#include "calc_instance.h"
#include "phylo_util.h"
#include <stdlib.h>
#include <string.h>
#include <libhmsbeagle/beagle.h>
#include "internal_like_calc_env.h"


char * convertBeagleEnumToCString(long beagleFlags, char * buffer) {
    if (buffer == 0L)
        return 0L;
    buffer[0] = '\0';
    if (beagleFlags & BEAGLE_FLAG_PRECISION_SINGLE) {
        if (buffer[0] != '\0')
            strcat(buffer, " | ");
        strcat(buffer, "BEAGLE_FLAG_PRECISION_SINGLE");
    }
    if (beagleFlags & BEAGLE_FLAG_PRECISION_DOUBLE) {
        if (buffer[0] != '\0')
            strcat(buffer, " | ");
        strcat(buffer, "BEAGLE_FLAG_PRECISION_DOUBLE");
    }
    if (beagleFlags & BEAGLE_FLAG_COMPUTATION_SYNCH) {
        if (buffer[0] != '\0')
            strcat(buffer, " | ");
        strcat(buffer, "BEAGLE_FLAG_COMPUTATION_SYNCH");
    }
    if (beagleFlags & BEAGLE_FLAG_COMPUTATION_ASYNCH) {
        if (buffer[0] != '\0')
            strcat(buffer, " | ");
        strcat(buffer, "BEAGLE_FLAG_COMPUTATION_ASYNCH");
    }
    if (beagleFlags & BEAGLE_FLAG_EIGEN_REAL) {
        if (buffer[0] != '\0')
            strcat(buffer, " | ");
        strcat(buffer, "BEAGLE_FLAG_EIGEN_REAL");
    }
    if (beagleFlags & BEAGLE_FLAG_EIGEN_COMPLEX) {
        if (buffer[0] != '\0')
            strcat(buffer, " | ");
        strcat(buffer, "BEAGLE_FLAG_EIGEN_COMPLEX");
    }
    if (beagleFlags & BEAGLE_FLAG_SCALING_MANUAL) {
        if (buffer[0] != '\0')
            strcat(buffer, " | ");
        strcat(buffer, "BEAGLE_FLAG_SCALING_MANUAL");
    }
    if (beagleFlags & BEAGLE_FLAG_SCALING_AUTO) {
        if (buffer[0] != '\0')
            strcat(buffer, " | ");
        strcat(buffer, "BEAGLE_FLAG_SCALING_AUTO");
    }
    if (beagleFlags & BEAGLE_FLAG_SCALING_ALWAYS) {
        if (buffer[0] != '\0')
            strcat(buffer, " | ");
        strcat(buffer, "BEAGLE_FLAG_SCALING_ALWAYS");
    }
    if (beagleFlags & BEAGLE_FLAG_SCALING_DYNAMIC) {
        if (buffer[0] != '\0')
            strcat(buffer, " | ");
        strcat(buffer, "BEAGLE_FLAG_SCALING_DYNAMIC");
    }
    if (beagleFlags & BEAGLE_FLAG_SCALERS_RAW) {
        if (buffer[0] != '\0')
            strcat(buffer, " | ");
        strcat(buffer, "BEAGLE_FLAG_SCALERS_RAW");
    }
    if (beagleFlags & BEAGLE_FLAG_SCALERS_LOG) {
        if (buffer[0] != '\0')
            strcat(buffer, " | ");
        strcat(buffer, "BEAGLE_FLAG_SCALERS_LOG");
    }
    if (beagleFlags & BEAGLE_FLAG_VECTOR_SSE) {
        if (buffer[0] != '\0')
            strcat(buffer, " | ");
        strcat(buffer, "BEAGLE_FLAG_VECTOR_SSE");
    }
    if (beagleFlags & BEAGLE_FLAG_VECTOR_NONE) {
        if (buffer[0] != '\0')
            strcat(buffer, " | ");
        strcat(buffer, "BEAGLE_FLAG_VECTOR_NONE");
    }
    if (beagleFlags & BEAGLE_FLAG_THREADING_OPENMP) {
        if (buffer[0] != '\0')
            strcat(buffer, " | ");
        strcat(buffer, "BEAGLE_FLAG_THREADING_OPENMP");
    }
    if (beagleFlags & BEAGLE_FLAG_THREADING_NONE) {
        if (buffer[0] != '\0')
            strcat(buffer, " | ");
        strcat(buffer, "BEAGLE_FLAG_THREADING_NONE");
    }
    if (beagleFlags & BEAGLE_FLAG_PROCESSOR_CPU) {
        if (buffer[0] != '\0')
            strcat(buffer, " | ");
        strcat(buffer, "BEAGLE_FLAG_PROCESSOR_CPU");
    }
    if (beagleFlags & BEAGLE_FLAG_PROCESSOR_GPU) {
        if (buffer[0] != '\0')
            strcat(buffer, " | ");
        strcat(buffer, "BEAGLE_FLAG_PROCESSOR_GPU");
    }
    if (beagleFlags & BEAGLE_FLAG_PROCESSOR_FPGA) {
        if (buffer[0] != '\0')
            strcat(buffer, " | ");
        strcat(buffer, "BEAGLE_FLAG_PROCESSOR_FPGA");
    }
    if (beagleFlags & BEAGLE_FLAG_PROCESSOR_CELL) {
        if (buffer[0] != '\0')
            strcat(buffer, " | ");
        strcat(buffer, "BEAGLE_FLAG_PROCESSOR_CELL");
    }
    return buffer;
}


BeagleOperation partialOperation(int destPartial, int outRescaler, int inRescaler, int leftPartial, int leftPrMat, int rightPartial, int rightPrMat) {
    BeagleOperation BOp;
    BOp.destinationPartials = destPartial;
	BOp.destinationScaleWrite = outRescaler;
	BOp.destinationScaleRead = inRescaler;
	BOp.child1Partials = leftPartial;
	BOp.child1TransitionMatrix = leftPrMat;
	BOp.child2Partials = rightPartial;
	BOp.child2TransitionMatrix = rightPrMat ; 
	return BOp;
}

#if !defined(BUILDING_FOR_PYTHON)
    void PyErr_NoMemory() {
        PYTBEAGLEHON_DEBUG_PRINTF("PyErr_NoMemory");
    }
    void PyErr_SetString(int c, const char *m) {
        PYTBEAGLEHON_DEBUG_PRINTF2("PyErr_SetString(%d, %s)", c, m);
    }
#endif

static struct LikeCalculatorInstance ** gAllInstances = 0;
static unsigned gLenAllInstancesArray = 0;

struct LikeCalculatorInstance * getLikeCalculatorInstance(long handle) {
	if (handle < gLenAllInstancesArray);
        return gAllInstances[handle];
    return 0L;
}


/** `numInstRateModels` is filled on output. */
const DSCTModelObj ** getModelList(long instanceHandle, unsigned int * numModels) {
    struct LikeCalculatorInstance * lci =  getLikeCalculatorInstance(instanceHandle);
    if (lci == 0L) {
        if (numModels)
            *numModels = 0;
        return 0L;
    }
    if (numModels)
        *numModels = lci->numInstRateModels;
    return (const DSCTModelObj **) lci->probModelArray;
}


long createNewLikeCalculatorInstance(void);
void freeLikeCalcInstanceFields(struct LikeCalculatorInstance *);
long allocateLikeCalcInstanceFields(struct LikeCalculatorInstance *, const ASRVObj ** asrvAliasForEachModel);






long createNewLikeCalculatorInstance(void) {
	/* In a multithreaded world, this should lock
		struct LikeCalculatorInstance ** gAllInstances = 0;
		unsigned gLenAllInstancesArray = 0;
	*/


	struct LikeCalculatorInstance ** tmp;
	struct LikeCalculatorInstance * newInstance;
	int i, destination;
	long newSize;
	if (gLenAllInstancesArray == 0) {
		tmp = (struct LikeCalculatorInstance **) mallocZeroedPointerArray(2);
		if (tmp == 0)
			goto errorExit;
		gLenAllInstancesArray = 2;
		gAllInstances = tmp;
	}

	/* Find an empty spot in the array */
	destination = -1;
	for (destination = 0; destination < gLenAllInstancesArray; ++destination) {
		if (gAllInstances[destination] == 0L)
			break;
	}
	if (destination >= gLenAllInstancesArray) {
		newSize = 2*gLenAllInstancesArray;
		tmp = (struct LikeCalculatorInstance **) mallocZeroedPointerArray(newSize);
		if (tmp == 0L)
			goto errorExit;
		for (i = 0; i < gLenAllInstancesArray; ++i)
			tmp[i] = gAllInstances[i];
		gAllInstances = tmp;
		destination = gLenAllInstancesArray;
		gLenAllInstancesArray = newSize;
	}
	assert (destination < gLenAllInstancesArray);
	assert (gAllInstances[destination] == 0L);

	newInstance = (struct LikeCalculatorInstance *)malloc(sizeof(struct LikeCalculatorInstance));
	PYTBEAGLEHON_DEBUG_PRINTF2("New struct LikeCalculatorInstance at %ld stored in global array %d\n", (long) newInstance, destination);
	if (newInstance == 0L)
		goto errorExit;
	gAllInstances[destination] = newInstance;
	return (long) destination;
	errorExit:
		PYTBEAGLEHON_DEBUG_PRINTF("In createNewLikeCalculatorInstance errorExit\n");
		return -1L;
}



void zeroLikeCalcInstanceFields(struct LikeCalculatorInstance * inst) {
    if (inst == 0L)
        return;
    inst->patternWeights = 0L;
	inst->probModelArray = 0L;
	inst->eigenSolutionStructs = 0L;
	inst->asrvAliasForEachModel = 0L;
	inst->beagleInstanceCreated = 0;
    inst->edgeLenScratch = 0L; 
    inst->probMatIndexScratch = 0L; 
    inst->probMatScratch = 0L;
    inst->stateCodeArrayScratch = 0L;
    inst->opScratch = 0L;
    inst->waitPartialIndexScratch = 0L;
    inst->categWeightIndexScratch = 0L;
    inst->categWeightScratch = 0L;
    inst->rootPartialIndexScratch = 0L;
    inst->stateFreqIndexScratch = 0L;
    inst->rootRescalerIndexScratch = 0L;

}

#if defined(BEAGLE_API_TRACE_PRINTING) && BEAGLE_API_TRACE_PRINTING
    int gTraceNumPatterns = 0;
#endif
/* returns 0 if memory allocation fails */
long allocateLikeCalcInstanceFields(struct LikeCalculatorInstance * t, const ASRVObj ** asrvAliasForEachModel) {
	unsigned int i;
	int doubleScratch, rc;
	int * resourceListPtr = 0L;
	int resourceListLen = 0;
	unsigned numRateCategoriesOnBeagle = 1;
	double rateOfOne = 1.0;
	if (!t) {
		PYTBEAGLEHON_DEBUG_PRINTF("EMPTY struct LikeCalculatorInstance in allocateLikeCalcInstanceFields\n");
		return 0;
	}
    t->patternWeights = 0L;
    if (t->numPatterns > 0) {
        t->patternWeights = (double *) malloc(t->numPatterns*sizeof(double));
        if (t->patternWeights == 0) {
            PYTBEAGLEHON_DEBUG_PRINTF("Could not alloc patternWeights in allocateLikeCalcInstanceFields\n");
		    goto errorExit;
        }
    }
	t->asrvAliasForEachModel = (const ASRVObj **) mallocZeroedPointerArray(t->numInstRateModels);
	if (t->asrvAliasForEachModel == 0){
		PYTBEAGLEHON_DEBUG_PRINTF("Could not alloc asrvAliasForEachModel in allocateLikeCalcInstanceFields\n");
		goto errorExit;
	}
	/* each model is a DSCTModelObj which is a PyObject */
	t->probModelArray = 0L;
	t->probModelArray = (DSCTModelObj **) mallocZeroedPointerArray(t->numInstRateModels);
	if (t->probModelArray == 0){
		PYTBEAGLEHON_DEBUG_PRINTF("Could not alloc probModelArray in allocateLikeCalcInstanceFields\n");
		goto errorExit;
	}
    t->numRateCategories = 0;
    for (i = 0; i < t->numInstRateModels; ++i) {
		t->probModelArray[i] = dsctModelNew(t->numStates);
		if (t->probModelArray[i] == 0L) {
			PYTBEAGLEHON_DEBUG_PRINTF1("Could not alloc probModelArray[%d] in allocateLikeCalcInstanceFields\n", i);
			goto errorExit;
		}
		Py_INCREF(t->probModelArray[i]);
		t->probModelArray[i]->likeCalcInstanceAlias = t;
        if (asrvAliasForEachModel != 0) {
            t->asrvAliasForEachModel[i] = asrvAliasForEachModel[i];
            if (t->asrvAliasForEachModel[i])
                t->numRateCategories += t->asrvAliasForEachModel[i]->n;
            else
                t->numRateCategories += 1; /* no rate het */
        }
        else {
            t->asrvAliasForEachModel[i] = 0L;
            t->numRateCategories += 1; /* no rate het */
        }
	}

	/* for each eigensolution, we need a EigenSolutionStruct
		currently we will also create a EigenCalcScratchpad for each EigenSolutionStruct
	*/
	t->eigenSolutionStructs = 0L;
	t->eigenSolutionStructs = (EigenSolutionStruct **) mallocZeroedPointerArray(t->numEigenStorage);
	if (t->eigenSolutionStructs == 0) {
		PYTBEAGLEHON_DEBUG_PRINTF("Could not alloc eigenSolutionStructs in allocateLikeCalcInstanceFields\n");
		goto errorExit;
	}
	for (i = 0; i < t->numEigenStorage; ++i) {
		t->eigenSolutionStructs[i] = eigenSolutionStructNew(t->numStates);
		if (t->eigenSolutionStructs[i] == 0L) {
			PYTBEAGLEHON_DEBUG_PRINTF1("Could not alloc eigenSolutionStructs[%d] in allocateLikeCalcInstanceFields\n", i);
			goto errorExit;
		}
		t->eigenSolutionStructs[i]->beagleEigenBufferIndex = i;
	}
    
    t->edgeLenScratch = (double *)malloc(t->numProbMats*sizeof(double));
	if (t->edgeLenScratch == 0) {
		PYTBEAGLEHON_DEBUG_PRINTF("Could not alloc edgeLenScratch in allocateLikeCalcInstanceFields\n");
		goto errorExit;
	}
    t->probMatIndexScratch = (int *)malloc(t->numProbMats*sizeof(int ));
	if (t->probMatIndexScratch == 0) {
		PYTBEAGLEHON_DEBUG_PRINTF("Could not alloc probMatIndexScratch in allocateLikeCalcInstanceFields\n");
		goto errorExit;
	}
    t->probMatScratch = allocateDblMatrix(t->numStates, t->numStates);
	if (t->probMatScratch == 0) {
		PYTBEAGLEHON_DEBUG_PRINTF("Could not alloc probMatScratch in allocateLikeCalcInstanceFields\n");
		goto errorExit;
	}
    t->stateCodeArrayScratch = (int *)malloc(t->numPatterns*sizeof(int));
	if (t->stateCodeArrayScratch == 0) {
		PYTBEAGLEHON_DEBUG_PRINTF("Could not alloc stateCodeArrayScratch in allocateLikeCalcInstanceFields\n");
		goto errorExit;
	}
    t->opScratch = (BeagleOperation *)malloc(t->numPartialStructs*sizeof(BeagleOperation));
	if (t->opScratch == 0) {
		PYTBEAGLEHON_DEBUG_PRINTF("Could not alloc opScratch in allocateLikeCalcInstanceFields\n");
		goto errorExit;
	}
    t->waitPartialIndexScratch = (int *)malloc(t->numPartialStructs*sizeof(int));
	if (t->waitPartialIndexScratch == 0) {
		PYTBEAGLEHON_DEBUG_PRINTF("Could not alloc waitPartialIndexScratch in allocateLikeCalcInstanceFields\n");
		goto errorExit;
	}
	
	/* this next two should not really be numEigenStorage long, but that is what beagle wants
	   we make it 4 times too long so that we can use the array for:
	    categWeightIndexScratch
	    rootPartialIndexScratch
	    stateFreqIndexScratch, and 
	    rootRescalerIndexScratch
    */
    t->categWeightIndexScratch = (int *)malloc(4*t->numEigenStorage*sizeof(int));
	if (t->categWeightIndexScratch == 0) {
		PYTBEAGLEHON_DEBUG_PRINTF("Could not alloc categWeightIndexScratch in allocateLikeCalcInstanceFields\n");
		goto errorExit;
	}
    t->rootPartialIndexScratch = t->categWeightIndexScratch + (t->numEigenStorage);
    t->stateFreqIndexScratch = t->categWeightIndexScratch + 2*(t->numEigenStorage);
    t->rootRescalerIndexScratch = t->categWeightIndexScratch + 3*(t->numEigenStorage);

	
	doubleScratch = (t->numEigenStorage > t->numStates ? t->numEigenStorage : t->numStates);
    t->categWeightScratch = (double *)malloc(doubleScratch*sizeof(double));
	if (t->categWeightScratch == 0) {
		PYTBEAGLEHON_DEBUG_PRINTF("Could not alloc categWeightScratch in allocateLikeCalcInstanceFields\n");
		goto errorExit;
	}

    BeagleResourceList * brl = beagleGetResourceList();
    if (brl == 0) {
        PYTBEAGLEHON_DEBUG_PRINTF("Could not alloc get beagle resource list in allocateLikeCalcInstanceFields\n");
        goto errorExit;
    }
    if (t->resourceIndex > brl->length) {
        PYTBEAGLEHON_DEBUG_PRINTF2("Resource %d is out of range (length = %d)\n", t->resourceIndex, brl->length);
        goto errorExit;
    }
    if (t->resourceIndex >= 0) {
        resourceListPtr = &(t->resourceIndex);
        resourceListLen = 1;
    }
    
    BeagleInstanceDetails beagleInstanceDetails;
#   if defined(BEAGLE_API_TRACE_PRINTING) && BEAGLE_API_TRACE_PRINTING
        PYTBEAGLEHON_DEBUG_PRINTF1("/* BEAGLE_API Call */ resourceIndex = %d;\n", t->resourceIndex);
        char prefStr[1000];
        char reqStr[1000];
        PYTBEAGLEHON_DEBUG_PRINTF2("/* BEAGLE_API Call */ resourcePref = %s ;\n/* BEAGLE_API Call */ resourceReq = %s ; \n", convertBeagleEnumToCString(t->resourcePref, prefStr), convertBeagleEnumToCString(t->resourceReq, reqStr));
        PYTBEAGLEHON_DEBUG_PRINTF4("/* BEAGLE_API Call */ rc = beagleCreateInstance(%d, %d, %d, %d", t->numLeaves, t->numPartialStructs, t->numStateCodeArrays, t->numStates);
        PYTBEAGLEHON_DEBUG_PRINTF4(", %ld, %d, %d, 1, %d, &resourceIndex, 1, resourcePref, resourceReq, &beagleInstanceDetails);\n", t->numPatterns, t->numEigenStorage, t->numProbMats, t->numRescalingsMultipliers);
#   endif
    
    rc = beagleCreateInstance(t->numLeaves,
                                                  t->numPartialStructs,
                                                  t->numStateCodeArrays,
                                                  t->numStates,
                                                  t->numPatterns,
                                                  t->numEigenStorage,
                                                  t->numProbMats,
                                                  1, /* we take care of the asrv at a higher level t->numRateCategories, */
                                                  t->numRescalingsMultipliers,
                                                  resourceListPtr,
                                                  resourceListLen,
                                                  t->resourcePref,
                                                  t->resourceReq,
                                                  &beagleInstanceDetails);
	PYTBEAGLEHON_DEBUG_PRINTF1("beagleCreateInstance returned %d\n", rc);
    if (rc != BEAGLE_SUCCESS) {
		PYTBEAGLEHON_DEBUG_PRINTF1("beagleCreateInstance failed with error code %d\n", rc);
		goto errorExit;
    }
    
#   if defined(BEAGLE_API_TRACE_PRINTING) && BEAGLE_API_TRACE_PRINTING
        PYTBEAGLEHON_DEBUG_PRINTF2("/* BEAGLE_API Call */ double * inst%dpatternWeights = (double *) malloc(%ld*sizeof(double));\n", rc, t->numPatterns);
        PYTBEAGLEHON_DEBUG_PRINTF2("/* BEAGLE_API Call */ double * inst%dedgeLenArray = (double *) malloc(%d*sizeof(double));\n", rc, t->numProbMats);
        PYTBEAGLEHON_DEBUG_PRINTF2("/* BEAGLE_API Call */ double * inst%dStateFreq = (double *) malloc(%d*sizeof(double));\n", rc, t->numEigenStorage);
        PYTBEAGLEHON_DEBUG_PRINTF2("/* BEAGLE_API Call */ double * inst%dCategWeight = (double *) malloc(%d*sizeof(double));\n", rc, t->numEigenStorage);
        PYTBEAGLEHON_DEBUG_PRINTF2("/* BEAGLE_API Call */ int * inst%dCategWeightIndex = (int *) malloc(%d*sizeof(int));\n", rc, t->numEigenStorage);
        PYTBEAGLEHON_DEBUG_PRINTF2("/* BEAGLE_API Call */ int * inst%dStateFreqIndex = (int *) malloc(%d*sizeof(int));\n", rc, t->numEigenStorage);
        PYTBEAGLEHON_DEBUG_PRINTF2("/* BEAGLE_API Call */ int * inst%dRootPartialIndex = (int *) malloc(%d*sizeof(int));\n", rc, t->numEigenStorage);
        PYTBEAGLEHON_DEBUG_PRINTF2("/* BEAGLE_API Call */ int * inst%dRootRescalerIndex = (int *) malloc(%d*sizeof(int));\n", rc, t->numEigenStorage);
        PYTBEAGLEHON_DEBUG_PRINTF2("/* BEAGLE_API Call */ int * inst%dpatternWeights = (int *) malloc(%d*sizeof(int));\n", rc, t->numProbMats);
        PYTBEAGLEHON_DEBUG_PRINTF2("/* BEAGLE_API Call */ int * inst%dWait = (int *) malloc(%d*sizeof(int));\n", rc, t->numPartialStructs);
        PYTBEAGLEHON_DEBUG_PRINTF3("/* BEAGLE_API Call */ BeagleOperation * inst%dOp = (BeagleOperation *) malloc(%d*sizeof(BeagleOperation)); double inst%dLnL;\n", rc, t->numPartialStructs, rc);
        PYTBEAGLEHON_DEBUG_PRINTF3("/* BEAGLE_API Call */ double ** inst%dEigenVec = allocateDblMatrix(%d, %d);\n", rc, t->numStates, t->numStates);
        PYTBEAGLEHON_DEBUG_PRINTF3("/* BEAGLE_API Call */ double ** inst%dInvEigenVec = allocateDblMatrix(%d, %d);\n", rc, t->numStates, t->numStates);
        PYTBEAGLEHON_DEBUG_PRINTF2("/* BEAGLE_API Call */ double * inst%dEigenVal = (double *) malloc(%d*sizeof(double));\n", rc, t->numStates);
        gTraceNumPatterns = t->numPatterns;
#   endif

    
    
    t->beagleInstanceIndex = rc;
	t->beagleInstanceCreated = 1;
    if (numRateCategoriesOnBeagle == 1) {
#       if defined(BEAGLE_API_TRACE_PRINTING) && BEAGLE_API_TRACE_PRINTING
            PYTBEAGLEHON_DEBUG_PRINTF3("/* BEAGLE_API Call */ double * inst%dRate = 1.0;  beagleSetCategoryRates(%d, &inst%dRate);\n", rc, rc, rc);
#       endif
        beagleSetCategoryRates(t->beagleInstanceIndex, &rateOfOne);
    }

	return 1L;
	errorExit:
		PYTBEAGLEHON_DEBUG_PRINTF("allocateLikeCalcInstanceFields errorExit\n");
		freeLikeCalcInstanceFields(t);
		return 0L;
}



int freeLikeCalculatorInstance(long handle) {
	struct LikeCalculatorInstance * inst = getLikeCalculatorInstance(handle);
	PYTBEAGLEHON_DEBUG_PRINTF2("freeLikeCalculatorInstance(%ld) -> %ld\n", handle, (long) inst);
	if (inst) {
		freeLikeCalcInstanceFields(inst);
		free(inst);
	}
	else
	    return BEAGLE_ERROR_OUT_OF_RANGE;
	gAllInstances[handle] = 0L;
	return 0;
}


/* public methods below here */

long createLikelihoodCalcInstance(
        unsigned int numLeaves,  /* I'm not sure why beagle needs to know this, but... */
        unsigned long numPatterns, /* the number of data patterns in this subset */
        const double * patternWeights, /* array of weights.  Length must be numPatterns */
        unsigned int numStates, /* the number of states in the data */
        unsigned int numStateCodeArrays, /* The number of data structs to allocate for storing leaf data. This is often equal to the number of leaves. If you have partial ambiguity in a tip, then you'll need to use a "partial" rather than a "state code" array for its storage.*/
        unsigned int numPartialStructs,  /* The number of data structures to allocate to store partial likelihood arrays. */
        unsigned int numInstRateModels, /* the number of distinct models (q-matrices) to allocated */
        const ASRVObj ** asrvAliasForEachModel, /* 0L or an array of length numInstRateModels. 0L indicates no rate heterogeneity for any model. Otherwise this should be an array of length numInstRateModels of pointers to the asrv object used for each model.  This will determine the number of probability matrices stored in the calculator instance */
        const unsigned int numProbMatsToAllocate, /* Number of transition probability matrices to allocate */ 
        unsigned int numEigenStorageStructs, /* number of eigen solution storage structures to allocate -- should be as large as the number of numInstRateModels (unless you are going to recalculate the eigensolution each time) */ 
        unsigned int numRescalingsMultipliers, /* the number of arrays to allocate for avoiding underflow */
        int resourceIndex, /* the index of the computational resource to use */
        long resourcePref,
        long resourceReq) {
	struct LikeCalculatorInstance * calcInstancePtr;
	long handle = createNewLikeCalculatorInstance();
	if (handle >= 0) {
		calcInstancePtr = getLikeCalculatorInstance(handle);
		assert(calcInstancePtr != 0L);
		zeroLikeCalcInstanceFields(calcInstancePtr);
		calcInstancePtr->numLeaves = numLeaves;
		calcInstancePtr->numPatterns = numPatterns;
		calcInstancePtr->numStates = numStates;
		calcInstancePtr->numStateCodeArrays = numStateCodeArrays;
		calcInstancePtr->numPartialStructs = numPartialStructs;
		calcInstancePtr->numInstRateModels = numInstRateModels;
		calcInstancePtr->numProbMats = numProbMatsToAllocate;
		calcInstancePtr->numEigenStorage = numEigenStorageStructs;
		calcInstancePtr->numRescalingsMultipliers = numRescalingsMultipliers;
		calcInstancePtr->resourceIndex = resourceIndex;
		calcInstancePtr->resourcePref = resourcePref;
		calcInstancePtr->resourceReq = resourceReq;
		calcInstancePtr->probModelArray = 0L;

		calcInstancePtr->eigenSolutionStructs = 0L;
		if (allocateLikeCalcInstanceFields(calcInstancePtr, asrvAliasForEachModel) == 0L) {
			PYTBEAGLEHON_DEBUG_PRINTF("Deleting calcInstancePtr that was just created because of allocateLikeCalcInstanceFields failure\n");
			freeLikeCalculatorInstance(handle);
			return BEAGLE_ERROR_OUT_OF_MEMORY;
		}
		
		setPatternWeights(handle, patternWeights);
	}
	return handle;
}


int setPatternWeights(long handle, const double * patternWeights) {
    unsigned i;
    struct LikeCalculatorInstance * calcInstancePtr = getLikeCalculatorInstance(handle);
    if (calcInstancePtr == 0L);
        return BEAGLE_ERROR_OUT_OF_RANGE;
    
    /* Set pattern weights and send them to beagle */
    for (i = 0; i < calcInstancePtr->numPatterns; ++i) {
        if (patternWeights == 0L)
            calcInstancePtr->patternWeights[i] = 1.0;
        else
            calcInstancePtr->patternWeights[i] = patternWeights[i];
    }
#   if defined(BEAGLE_API_TRACE_PRINTING) && BEAGLE_API_TRACE_PRINTING
        PYTBEAGLEHON_DEBUG_PRINTF("/* BEAGLE_API Call */ ");
        for (i = 0; i < calcInstancePtr->numPatterns; ++i) {
            PYTBEAGLEHON_DEBUG_PRINTF3("inst%dpatternWeights[%u] = %f; ", calcInstancePtr->beagleInstanceIndex, i,  calcInstancePtr->patternWeights[i]);
        }
        PYTBEAGLEHON_DEBUG_PRINTF("\n");
#   endif
    beagleSetPatternWeights(calcInstancePtr->beagleInstanceIndex, calcInstancePtr->patternWeights);
    return 0;
}


void freeLikeCalcInstanceFields(struct LikeCalculatorInstance * inst) {
	unsigned int i;
	if (inst == 0L)
		return;
    if (inst->beagleInstanceCreated) {
#       if defined(BEAGLE_API_TRACE_PRINTING) && BEAGLE_API_TRACE_PRINTING
            PYTBEAGLEHON_DEBUG_PRINTF1("/* BEAGLE_API Call */ beagleFinalizeInstance(%d);\n", inst->beagleInstanceIndex);
#       endif
        beagleFinalizeInstance(inst->beagleInstanceIndex);
    }

	if (inst->patternWeights != 0L)
	    free(inst->patternWeights);
    inst->patternWeights = 0L;
    
	for (i = 0; i < inst->numInstRateModels; ++i) {
		if (inst->probModelArray[i] != 0L) {
			Py_DECREF(inst->probModelArray[i]);
			/* cdsctm_dtor(inst->probModelArray[i]); */
			inst->probModelArray[i] = 0L;
		}
#       if defined (BUILDING_FOR_PYTHON)
        if (inst->asrvAliasForEachModel && inst->asrvAliasForEachModel[i]) {
            Py_DECREF((PyObject *) inst->asrvAliasForEachModel[i]);
        }
#       endif

	}
	PYTBEAGLEHON_DEBUG_PRINTF("freeing struct LikeCalculatorInstance->probModelArray...\n");
	if (inst->probModelArray)
	    free(inst->probModelArray);
	inst->probModelArray = 0L;

	PYTBEAGLEHON_DEBUG_PRINTF("freeing struct LikeCalculatorInstance->eigenSolutionStructs elements...\n");
	for (i = 0; i < inst->numEigenStorage; ++i) {
		eigenSolutionStructDtor(inst->eigenSolutionStructs[i]);
		inst->eigenSolutionStructs[i] = 0L;
	}

	if (inst->eigenSolutionStructs)
	    free(inst->eigenSolutionStructs);
	inst->eigenSolutionStructs = 0L;

	inst->numEigenStorage = 0;
	
	if (inst->probMatIndexScratch)
	    free(inst->probMatIndexScratch);
	inst->probMatIndexScratch = 0L;

	if (inst->edgeLenScratch)
	    free(inst->edgeLenScratch);
	inst->edgeLenScratch = 0L;

	freeDblMatrix(inst->probMatScratch);
	inst->probMatScratch = 0L;

    if (inst->stateCodeArrayScratch)
	    free(inst->stateCodeArrayScratch);
	inst->stateCodeArrayScratch = 0L;

	if (inst->opScratch)
	    free(inst->opScratch);
	inst->opScratch = 0L;

	if (inst->waitPartialIndexScratch)
	    free(inst->waitPartialIndexScratch);
	inst->waitPartialIndexScratch = 0L;

	if (inst->categWeightIndexScratch)
	    free(inst->categWeightIndexScratch);
	inst->categWeightIndexScratch = 0L;
    inst->rootPartialIndexScratch = 0L;
    inst->stateFreqIndexScratch = 0L;
    inst->rootRescalerIndexScratch = 0L;

	if (inst->categWeightScratch)
	    free(inst->categWeightScratch);
	inst->categWeightScratch = 0L;

}



unsigned getNumComputationalResources() {
#   if defined(BEAGLE_API_TRACE_PRINTING) && BEAGLE_API_TRACE_PRINTING
        PYTBEAGLEHON_DEBUG_PRINTF("/* BEAGLE_API Call */ beagleGetResourceList();\n");
#   endif
    BeagleResourceList *brl =  beagleGetResourceList();
    if (brl)
        return (unsigned) brl->length;
    return 0;
}


int getComputationalResourceDetails(int resourceIndex,
                                          char * resourceName, /* must be an array of at least 81 chars. This will be filled in with the first 80 characters of the resource name */
                                          char * description, /* must be an array of at least 81 chars. This will be filled in with the first 80 characters of the resource name */
                                          long * supportedFlags, /* Supported features (as beagle flags) */
                                          long * requiredFlags /* as beagle flags */
                                          )
    {
    const BeagleResourceList *brl =  beagleGetResourceList();
    if (brl == 0L)
        return BEAGLE_ERROR_GENERAL;
    if (resourceIndex >= (int) brl->length || resourceIndex < 0)
        return BEAGLE_ERROR_OUT_OF_RANGE;
    BeagleResource * br = brl->list + resourceIndex;
    if (br == 0L)
        return BEAGLE_ERROR_UNINITIALIZED_INSTANCE;


	if (! (br->supportFlags & BEAGLE_FLAG_SCALING_ALWAYS) ) /* \TEMP only supporting ALWAYS scaling */
	    return BEAGLE_ERROR_GENERAL;
	



    const unsigned maxIndex = 80;
    if (strlen(br->name) > maxIndex)
        br->name[maxIndex] = '\0';
    strcpy(resourceName, br->name);

    if (strlen(br->description) > maxIndex)
        br->description[maxIndex] = '\0';
    strcpy(description, br->description);
    
    if (supportedFlags) {  
        *supportedFlags = br->supportFlags;
        /* \TEMP only supporting ALWAYS scaling */
        *supportedFlags &= (~BEAGLE_FLAG_SCALING_MANUAL);
        *supportedFlags &= (~BEAGLE_FLAG_SCALING_AUTO);
        *supportedFlags &= (~BEAGLE_FLAG_SCALING_DYNAMIC);
        *supportedFlags |= BEAGLE_FLAG_SCALING_ALWAYS;
    }
    if (requiredFlags) {
        *requiredFlags = br->requiredFlags;
        /* \TEMP only supporting ALWAYS scaling */
        *requiredFlags &= (~BEAGLE_FLAG_SCALING_MANUAL);
        *requiredFlags &= (~BEAGLE_FLAG_SCALING_AUTO);
        *requiredFlags &= (~BEAGLE_FLAG_SCALING_DYNAMIC);
        *requiredFlags |= BEAGLE_FLAG_SCALING_DYNAMIC;
    }

    return 0;
    
}


int calcPrMats(long handle, 
               int eigenIndex,
               unsigned numToCalc,
               const double * edgeLenArray,
               const int * probMatIndexArray) {
    struct LikeCalculatorInstance * lci;
    EigenSolutionStruct * ess;
    int bEigen;
    lci = getLikeCalculatorInstance(handle);
    if (lci == 0L || numToCalc > lci->numProbMats) {
		return BEAGLE_ERROR_OUT_OF_RANGE;
    }
    if (eigenIndex < 0 || eigenIndex >= lci->numEigenStorage) {
		return BEAGLE_ERROR_OUT_OF_RANGE;
    }
    ess = lci->eigenSolutionStructs[eigenIndex];
    if (ess == 0L || ess->beagleEigenBufferIndex < 0) {
		return BEAGLE_ERROR_OUT_OF_RANGE;
    }
    bEigen = ess->beagleEigenBufferIndex;
    PYTBEAGLEHON_DEBUG_PRINTF5("In C. calling beagleUpdateTransitionMatrices(%d, %d, [%d...], 0, 0, [%lf...],%d)\n", (int) lci->beagleInstanceIndex, bEigen, *probMatIndexArray, *edgeLenArray, numToCalc); 
#   if defined(BEAGLE_API_TRACE_PRINTING) && BEAGLE_API_TRACE_PRINTING
        int i;
        PYTBEAGLEHON_DEBUG_PRINTF("/* BEAGLE_API Call */ ");
        for (i = 0; i < numToCalc; ++i) {
            PYTBEAGLEHON_DEBUG_PRINTF3("inst%dProbMatArray[%d] = %d; ", lci->beagleInstanceIndex, i, probMatIndexArray[i]);
            PYTBEAGLEHON_DEBUG_PRINTF3("inst%dedgeLenArray[%d] = %lf; ", lci->beagleInstanceIndex, i, edgeLenArray[i]);
        }
        PYTBEAGLEHON_DEBUG_PRINTF5("\n/* BEAGLE_API Call */ beagleUpdateTransitionMatrices(%d, %d, inst%dProbMatArray, inst%dedgeLenArray, %d);\n", lci->beagleInstanceIndex, bEigen, lci->beagleInstanceIndex, lci->beagleInstanceIndex, numToCalc);
#   endif
    return beagleUpdateTransitionMatrices(lci->beagleInstanceIndex,
                                   bEigen,
                                   probMatIndexArray,
                                   0L, /*firstDerivativeIndices*/
                                   0L, /*secondDerivativeIndices*/
                                   edgeLenArray,
                                   numToCalc);
}


int fetchPrMat(long handle, int probMatIndex, double * flattenedMat) {
    struct LikeCalculatorInstance * lci;
    int rc;
    lci = getLikeCalculatorInstance(handle);
    if (lci == 0L || probMatIndex >= lci->numProbMats) {
		return BEAGLE_ERROR_OUT_OF_RANGE;
    }
    PYTBEAGLEHON_DEBUG_PRINTF2("In C. calling beagleGetTransitionMatrix(%d, %d, ...)\n", (int) lci->beagleInstanceIndex, probMatIndex); 
#   if defined(BEAGLE_API_TRACE_PRINTING) && BEAGLE_API_TRACE_PRINTING
        PYTBEAGLEHON_DEBUG_PRINTF3("/* BEAGLE_API Call */ beagleGetTransitionMatrix(%d, %d, inst%dflattened)", lci->beagleInstanceIndex, probMatIndex, lci->beagleInstanceIndex);
#   endif
    rc = beagleGetTransitionMatrix(lci->beagleInstanceIndex,
								     probMatIndex,
								     flattenedMat);
    PYTBEAGLEHON_DEBUG_PRINTF2("In C. First two elements from beagleGetTransitionMatrix = %lf, %lf, ...\n" , flattenedMat[0], flattenedMat[1]); 
    
    return rc;
}

int setStateCodeArray(long handle, int stateCodeArrayIndex, const int * stateCodes) {
    struct LikeCalculatorInstance * lci;
    int rc;
    lci = getLikeCalculatorInstance(handle);
    if (lci == 0L || stateCodeArrayIndex >= lci->numStateCodeArrays) {
		return BEAGLE_ERROR_OUT_OF_RANGE;
    }
#   if defined(BEAGLE_API_TRACE_PRINTING) && BEAGLE_API_TRACE_PRINTING
        PYTBEAGLEHON_DEBUG_PRINTF("/* BEAGLE_API Call */ ");
        for (rc = 0; rc < lci->numPatterns; ++rc) {
            PYTBEAGLEHON_DEBUG_PRINTF3("inst%dstateCodeScratch[%d] = %d; ", lci->beagleInstanceIndex, rc, stateCodes[rc]);
        }
        PYTBEAGLEHON_DEBUG_PRINTF3("\n/* BEAGLE_API Call */ beagleGetTransitionMatrix(%d, %d, inst%dstateCodeScratch)", lci->beagleInstanceIndex, stateCodeArrayIndex, lci->beagleInstanceIndex);
#   endif
    rc = beagleSetTipStates(lci->beagleInstanceIndex, stateCodeArrayIndex, stateCodes);
    return rc;
}



int calcPartials(long handle, const BeagleOperation * opArray, unsigned numOps, const int * waitPartialIndex, int numPartialsToWaitFor) {
    struct LikeCalculatorInstance * lci;
    int rc = BEAGLE_SUCCESS;
    lci = getLikeCalculatorInstance(handle);
    if (lci == 0L || numOps >= lci->numPartialStructs) {
		return BEAGLE_ERROR_OUT_OF_RANGE;
    }
    if (numOps > 0) {
#       if defined(BEAGLE_API_TRACE_PRINTING) && BEAGLE_API_TRACE_PRINTING
            for (rc = 0; rc < numOps; ++rc) {
                PYTBEAGLEHON_DEBUG_PRINTF3("/* BEAGLE_API Call */ inst%dOp[%d].destinationPartials = %d; ", lci->beagleInstanceIndex, rc, opArray[rc].destinationPartials);
                PYTBEAGLEHON_DEBUG_PRINTF3("inst%dOp[%d].destinationScaleWrite = %d; ", lci->beagleInstanceIndex, rc, opArray[rc].destinationScaleWrite);
                PYTBEAGLEHON_DEBUG_PRINTF3("inst%dOp[%d].destinationScaleRead = %d; ", lci->beagleInstanceIndex, rc, opArray[rc].destinationScaleRead);
                PYTBEAGLEHON_DEBUG_PRINTF3("inst%dOp[%d].child1Partials = %d; ", lci->beagleInstanceIndex, rc, opArray[rc].child1Partials);
                PYTBEAGLEHON_DEBUG_PRINTF3("inst%dOp[%d].child1TransitionMatrix = %d; ", lci->beagleInstanceIndex, rc, opArray[rc].child1TransitionMatrix);
                PYTBEAGLEHON_DEBUG_PRINTF3("inst%dOp[%d].child2Partials = %d; ", lci->beagleInstanceIndex, rc, opArray[rc].child2Partials);
                PYTBEAGLEHON_DEBUG_PRINTF3("inst%dOp[%d].child2TransitionMatrix = %d;\n", lci->beagleInstanceIndex, rc, opArray[rc].child2TransitionMatrix);
            }
            PYTBEAGLEHON_DEBUG_PRINTF3("\n/* BEAGLE_API Call */ beagleUpdatePartials(%d, inst%dOp, %d, BEAGLE_OP_NONE);\n", lci->beagleInstanceIndex, lci->beagleInstanceIndex, numOps);
#       endif
        rc = beagleUpdatePartials(lci->beagleInstanceIndex, opArray, numOps, BEAGLE_OP_NONE); /*TEMP: I think BEAGLE_OP_NONE is appropriate if we are in always scale mode*/
        if (rc != BEAGLE_SUCCESS) {
            PYTBEAGLEHON_DEBUG_PRINTF("Error in beagleUpdatePartials");
            return rc;
        }
    }
    if (numPartialsToWaitFor >= 0) {
#       if defined(BEAGLE_API_TRACE_PRINTING) && BEAGLE_API_TRACE_PRINTING
            for (rc = 0; rc < numPartialsToWaitFor; ++rc) {
                PYTBEAGLEHON_DEBUG_PRINTF3("/* BEAGLE_API Call */ inst%dWait[%d] = %d;\n", lci->beagleInstanceIndex, rc, waitPartialIndex[rc]);
            }
            PYTBEAGLEHON_DEBUG_PRINTF3("\n/* BEAGLE_API Call */ beagleWaitForPartials(%d, inst%dWait, %d);\n", lci->beagleInstanceIndex, lci->beagleInstanceIndex, numPartialsToWaitFor);
#       endif
        rc = beagleWaitForPartials(lci->beagleInstanceIndex, waitPartialIndex, numPartialsToWaitFor);
        if (rc != BEAGLE_SUCCESS) {
            PYTBEAGLEHON_DEBUG_PRINTF("Error in beagleWaitForPartials");
            return rc;
        }
    }
    return rc;
}

int setSingletonCategoryWeights(long handle, const int * indexList, const double *wtList, int numCateg) {
    int i;
    struct LikeCalculatorInstance * lci;
    int rc = BEAGLE_SUCCESS;
    lci = getLikeCalculatorInstance(handle);
    if (lci == 0L || numCateg >= lci->numEigenStorage) {
		return BEAGLE_ERROR_OUT_OF_RANGE;
    }
    for (i = 0; i < numCateg; ++i) {
#       if defined(BEAGLE_API_TRACE_PRINTING) && BEAGLE_API_TRACE_PRINTING
            for (rc = 0; rc < numCateg; ++rc) {
                PYTBEAGLEHON_DEBUG_PRINTF3("/* BEAGLE_API Call */ inst%dCategWeightIndex[%d] = %d; ", lci->beagleInstanceIndex, rc, indexList[rc]);
                PYTBEAGLEHON_DEBUG_PRINTF3("inst%dCategWeightIndex[%d] = %lf;\n", lci->beagleInstanceIndex, rc, wtList[rc]);
            }
            PYTBEAGLEHON_DEBUG_PRINTF3("\n/* BEAGLE_API Call */ beagleSetCategoryWeights(%d, inst%dCategWeightIndex, inst%dCategWeight);\n", lci->beagleInstanceIndex, lci->beagleInstanceIndex, lci->beagleInstanceIndex);
#       endif
        rc = beagleSetCategoryWeights(lci->beagleInstanceIndex, indexList[i], &wtList[i]);
        if (rc != BEAGLE_SUCCESS) {
            PYTBEAGLEHON_DEBUG_PRINTF("Error in beagleSetCategoryWeights");
            return rc;
        }
    }
    return rc;
}



int setStateFreq(long handle, int bufferIndex, const double *freq) {
    struct LikeCalculatorInstance * lci;
    int rc = BEAGLE_SUCCESS;
    lci = getLikeCalculatorInstance(handle);
    if (lci == 0L || bufferIndex >= lci->numEigenStorage) {
		return BEAGLE_ERROR_OUT_OF_RANGE;
    }
#   if defined(BEAGLE_API_TRACE_PRINTING) && BEAGLE_API_TRACE_PRINTING
        for (rc = 0; rc < lci->numStates; ++rc) {
            PYTBEAGLEHON_DEBUG_PRINTF3("/* BEAGLE_API Call */ inst%dStateFreq[%d] = %lf;\n", lci->beagleInstanceIndex, rc, freq[rc]);
        }
        PYTBEAGLEHON_DEBUG_PRINTF3("/* BEAGLE_API Call */ beagleSetStateFrequencies(%d, %d, inst%dStateFreq);\n", lci->beagleInstanceIndex, bufferIndex, lci->beagleInstanceIndex);
#   endif
    rc = beagleSetStateFrequencies(lci->beagleInstanceIndex, bufferIndex, freq);
    if (rc != BEAGLE_SUCCESS) {
        PYTBEAGLEHON_DEBUG_PRINTF("Error in beagleSetStateFrequencies");
    }
    return rc;
}


int calcRootLnL(long handle, 
                const int * rootPartialIndex,
                const int * categWeightIndex,
                const int * stateFreqIndex, 
                const int * rootRescalerIndex, 
                int arrayLength,
                double * lnL) {
    struct LikeCalculatorInstance * lci;
    int rc = BEAGLE_SUCCESS;
    lci = getLikeCalculatorInstance(handle);
    if (lci == 0L || arrayLength >= lci->numEigenStorage) {
		return BEAGLE_ERROR_OUT_OF_RANGE;
    }
#   if defined(BEAGLE_API_TRACE_PRINTING) && BEAGLE_API_TRACE_PRINTING
        assert(arrayLength <= lci->numEigenStorage);
        for (rc = 0; rc < arrayLength; ++rc) {
            PYTBEAGLEHON_DEBUG_PRINTF3("/* BEAGLE_API Call */ inst%dRootPartialIndex[%d] = %d;\n", lci->beagleInstanceIndex, rc, rootPartialIndex[rc]);
            PYTBEAGLEHON_DEBUG_PRINTF3("/* BEAGLE_API Call */ inst%dCategWeightIndex[%d] = %d;\n", lci->beagleInstanceIndex, rc, categWeightIndex[rc]);
            PYTBEAGLEHON_DEBUG_PRINTF3("/* BEAGLE_API Call */ inst%dStateFreqIndex[%d] = %d;\n", lci->beagleInstanceIndex, rc, stateFreqIndex[rc]);
            PYTBEAGLEHON_DEBUG_PRINTF3("/* BEAGLE_API Call */ inst%dRootRescalerIndex[%d] = %d;\n", lci->beagleInstanceIndex, rc, rootRescalerIndex[rc]);
        }
        PYTBEAGLEHON_DEBUG_PRINTF4("/* BEAGLE_API Call */ beagleCalculateRootLogLikelihoods(%d, inst%dRootPartialIndex, inst%dCategWeightIndex, inst%dStateFreqIndex, ", lci->beagleInstanceIndex, lci->beagleInstanceIndex, lci->beagleInstanceIndex, lci->beagleInstanceIndex);
        PYTBEAGLEHON_DEBUG_PRINTF3("inst%dRootRescalerIndex, %d, &inst%dLnL);\n", lci->beagleInstanceIndex, lci->beagleInstanceIndex, lci->beagleInstanceIndex);
#   endif
    rc = beagleCalculateRootLogLikelihoods(lci->beagleInstanceIndex,
                                      rootPartialIndex,
                                      categWeightIndex,
                                      stateFreqIndex,
                                      rootRescalerIndex,
                                      arrayLength,
                                      lnL);
    if (rc != BEAGLE_SUCCESS) {
        PYTBEAGLEHON_DEBUG_PRINTF("Error in beagleCalculateRootLogLikelihoods");
    }
    return rc;
}




#   if defined(API_TRACE_PRINTING) && API_TRACE_PRINTING
        int getTraceModeModelIndex(struct LikeCalculatorInstance * LCI, DSCTModelObj * m) {
            int i;
            if (LCI == 0L || LCI->probModelArray == 0L)
                return -1;
            for (i = 0; i < LCI->numInstRateModels; ++ i) {
                if (LCI->probModelArray[i] == m)
                    return i;
            }
            return -1;
        }


#   endif


/**
  pytbeaglehon phylogenetic likelihood caluclations using beaglelib.

  Copyright 2010 Mark T. Holder (except where noted). All rights reserved.

  See "LICENSE.txt" for terms and conditions of usage.
*/
