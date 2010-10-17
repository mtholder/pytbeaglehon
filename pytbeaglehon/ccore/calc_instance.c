#include "calc_instance.h"
#include "phylo_util.h"
#include <stdlib.h>
#include <string.h>
#include <libhmsbeagle/beagle.h>
#include "internal_like_calc_env.h"


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

INLINE struct LikeCalculatorInstance * getLikeCalculatorInstance(long handle) {
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
}

/* returns 0 if memory allocation fails */
long allocateLikeCalcInstanceFields(struct LikeCalculatorInstance * t, const ASRVObj ** asrvAliasForEachModel) {
	unsigned int i;
	int * resourceListPtr = 0L;
	int resourceListLen = 0;
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
    t->beagleInstanceIndex = beagleCreateInstance(t->numLeaves,
                                                  t->numPartialStructs,
                                                  t->numStateCodeArrays,
                                                  t->numStates,
                                                  t->numPatterns,
                                                  t->numEigenStorage,
                                                  t->numProbMats,
                                                  t->numRateCategories,
                                                  t->numRescalingsMultipliers,
                                                  resourceListPtr,
                                                  resourceListLen,
                                                  t->resourcePref,
                                                  t->resourceReq,
                                                  &beagleInstanceDetails);
	t->beagleInstanceCreated = 1;
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
    beagleSetPatternWeights(calcInstancePtr->beagleInstanceIndex, calcInstancePtr->patternWeights);
    return 0;
}


void freeLikeCalcInstanceFields(struct LikeCalculatorInstance * inst) {
	if (inst == 0L)
		return;
    if (inst->beagleInstanceCreated)
        beagleFinalizeInstance(inst->beagleInstanceIndex);
	unsigned int i;
	if (inst->patternWeights != 0L)
	    free(inst->patternWeights);

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
	free(inst->probModelArray);
	inst->probModelArray = 0L;

	PYTBEAGLEHON_DEBUG_PRINTF("freeing struct LikeCalculatorInstance->eigenSolutionStructs elements...\n");
	for (i = 0; i < inst->numEigenStorage; ++i) {
		eigenSolutionStructDtor(inst->eigenSolutionStructs[i]);
		inst->eigenSolutionStructs[i] = 0L;
	}
	free(inst->eigenSolutionStructs);
	inst->eigenSolutionStructs = 0L;
	inst->numEigenStorage = 0;
}



unsigned getNumComputationalResources() {
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

    const unsigned maxIndex = 80;
    if (strlen(br->name) > maxIndex)
        br->name[maxIndex] = '\0';
    strcpy(resourceName, br->name);

    if (strlen(br->description) > maxIndex)
        br->description[maxIndex] = '\0';
    strcpy(description, br->description);
    
    if (supportedFlags)
        *supportedFlags = br->supportFlags;
    if (requiredFlags)
        *requiredFlags = br->requiredFlags;
    return 0;
    
}



/**
  pytbeaglehon phylogenetic likelihood caluclations using beaglelib.

  Copyright 2010 Mark T. Holder (except where noted). All rights reserved.

  See "LICENSE.txt" for terms and conditions of usage.
*/
