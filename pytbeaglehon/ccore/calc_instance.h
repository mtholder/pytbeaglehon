/**  Provides a handle for the memory allocated for likelihood calculations (other
    than the ASRV structurs
 */
#if ! defined(CALC_INSTANCE_H)
#define CALC_INSTANCE_H
#ifdef __cplusplus
extern "C" 
{
#endif



#include "pytbeaglehon_defs.h"
#include "asrv.h"


unsigned getNumComputationalResources(void);

/* returns 0 for success and a negtive number (a facet of BeagleReturnCodes for failure). */
int getComputationalResourceDetails(int resourceIndex,
                                          char * resourceName, /* must be an array of at least 81 chars. This will be filled in with the first 80 characters of the resource name */
                                          char * description, /* must be an array of at least 81 chars. This will be filled in with the first 80 characters of the resource name */
                                          long * supportedFlags, /* Supported features (as beagle flags) */
                                          long * requiredFlags /* as beagle flags */
                                          );

/* Beagle Flags
    BEAGLE_FLAG_PRECISION_SINGLE    = 1 << 0,    < Single precision computation 
    BEAGLE_FLAG_PRECISION_DOUBLE    = 1 << 1,    < Double precision computation 

    BEAGLE_FLAG_COMPUTATION_SYNCH   = 1 << 2,    < Synchronous computation (blocking) 
    BEAGLE_FLAG_COMPUTATION_ASYNCH  = 1 << 3,    < Asynchronous computation (non-blocking) 
    
    BEAGLE_FLAG_EIGEN_REAL          = 1 << 4,    < Real eigenvalue computation 
    BEAGLE_FLAG_EIGEN_COMPLEX       = 1 << 5,    < Complex eigenvalue computation 

    BEAGLE_FLAG_SCALING_MANUAL      = 1 << 6,    < Manual scaling 
    BEAGLE_FLAG_SCALING_AUTO        = 1 << 7,    < Auto-scaling on 
    BEAGLE_FLAG_SCALING_ALWAYS      = 1 << 8,    < Scale at every updatePartials 
    BEAGLE_FLAG_SCALING_DYNAMIC     = 1 << 19,   < Manual scaling with dynamic checking  
    
    BEAGLE_FLAG_SCALERS_RAW         = 1 << 9,    < Save raw scalers 
    BEAGLE_FLAG_SCALERS_LOG         = 1 << 10,   < Save log scalers 
    
    BEAGLE_FLAG_VECTOR_SSE          = 1 << 11,   < SSE computation 
    BEAGLE_FLAG_VECTOR_NONE         = 1 << 12,   < No vector computation 
    
    BEAGLE_FLAG_THREADING_OPENMP    = 1 << 13,   < OpenMP threading 
    BEAGLE_FLAG_THREADING_NONE      = 1 << 14,   < No threading 
    
    BEAGLE_FLAG_PROCESSOR_CPU       = 1 << 15,   < Use CPU as main processor 
    BEAGLE_FLAG_PROCESSOR_GPU       = 1 << 16,   < Use GPU as main processor 
    BEAGLE_FLAG_PROCESSOR_FPGA      = 1 << 17,   < Use FPGA as main processor 
    BEAGLE_FLAG_PROCESSOR_CELL      = 1 << 18    < Use Cell as main processor 
*/
/* Allocates a likelihood calculator of the appropriate size, and returns its handle
    negative return values indicate errors.
*/
long createLikelihoodCalcInstance(
        unsigned int numLeaves,  /* I'm not sure why beagle needs to know this, but... */
        unsigned long numPatterns, /* the number of data patterns in this subset */
        const double * patternWeights, /* array of weights.  Length must be numPatterns */
        unsigned int numStates, /* the number of states in the data */
        unsigned int numStateCodeArrays, /* The number of data structs to allocate for storing leaf data. This is often equal to the number of leaves. If you have partial ambiguity in a tip, then you'll need to use a "partial" rather than a "state code" array for its storage.*/
        unsigned int numPartialStructs,  /* The number of data structures to allocate to store partial likelihood arrays. */
        unsigned int numInstRateModels, /* the number of distinct models (q-matrices) to allocated */
        const ASRVObj ** asrvForEachModel, /* 0L or an array of length numInstRateModels. 0L indicates no rate heterogeneity for any model. Otherwise this should be an array of length numProbModels of pointers to the asrv object used for each model.  This will determine the number of probability matrices stored in the calculator instance */
        const unsigned int numProbMatsToAllocate, /* Number of transition probability matrices to allocate */ 
        unsigned int numEigenStorageStructs, /* number of eigen solution storage structures to allocate -- should be as large as the number of numInstRateModels (unless you are going to recalculate the eigensolution each time) */ 
        unsigned int numRescalingsMultipliers, /* the number of arrays to allocate for avoiding underflow */
        int resourceIndex, /* the index of the computational resource to use */
        long resourcePref, /* the beagle flags (see above) preferred in a computational resource */
        long resourceReq); /* the beagle flags (see above) required in a computational resource */

/* frees the likelihood calculator with the identifier `handle` 
    \returns 0 for success, and a BeagleReturnCode for an error
*/
int freeLikeCalculatorInstance(long handle);


/*  likeCalcHandle must be a handle obtained via createLikelihoodCalcInstance
\returns 0 or BeagleReturnCode for failure */
int setPatternWeights(long likeCalcHandle, const double * patternWeights);

/*  calculates `numToCalc` probability matrices from the eigen system stored 
    at `eigenIndex` using the branch lengths in the array `edgeLenArray` and stores
    the matrices in `probMatIndexArray`
\returns 0 or BeagleReturnCode for failure */
int calcPrMats(long handle, int eigenIndex, unsigned numToCalc, const double * edgeLenArray, const int * probMatIndexArray);
int fetchPrMat(long handle, int probMatIndex, double * flattenedMat);

int setStateCodeArray(long handle, int stateCodeArrayIndex, const int * state_codes);


/* needed for BeagleOperation type in calcPartials call */
#include <libhmsbeagle/beagle.h>

/** queues the partial operations in opArray.
    `waitPartialIndex` should be the index of a partial buffer to wait on (or < 0 to return immediatly) 
    \returns 0 or BeagleReturnCode for failure
*/
int calcPartials(long handle, const BeagleOperation * opArray, unsigned numOps, const int * waitPartialIndexList, int numPartialsToWaitFor);

int setSingletonCategoryWeights(long handle, const int * indexList, const double *wtList, int numCateg);
int setStateFreq(long handle, int bufferIndex, const double *freq);

int calcRootLnL(long handle, 
                const int * rootPartialIndex,
                const int * categWeightIndex,
                const int * stateFreqIndex, 
                const int * rootRescalerIndex, 
                int arrayLength,
                double * lnL);

#ifdef __cplusplus
}
/* end of extern c bit */
#endif






#endif

/**
  pytbeaglehon phylogenetic likelihood caluclations using beaglelib.

  Copyright 2010 Mark T. Holder (except where noted). All rights reserved.

  See "LICENSE.txt" for terms and conditions of usage.
*/
