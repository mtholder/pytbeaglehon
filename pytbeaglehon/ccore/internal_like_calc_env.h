/**  Provides a handle for the memory allocated for likelihood calculations (other
    than the ASRV structurs
 */
#if ! defined(INTERNAL_LIKE_CALC_ENV_H)
#define INTERNAL_LIKE_CALC_ENV_H

#include <libhmsbeagle/beagle.h>


#ifdef __cplusplus
extern "C" 
{
#endif

#include "discrete_state_model.h"


/*! Each LikeCalculatorInstance contains the internal storage needed to
 *      calculate likelihoods for one subset in a (possibly) partitioned matrix.
 *  Multiple models may be associated with the instance, but all model share the
 *      same number of states.
 *  All of the fields should be viewed as "read-only"
 */
struct LikeCalculatorInstance {
    int beagleInstanceIndex;
    int handleForSelf; /**< The handle that applies to this LikeCalculatorInstance (the arg to #getLikeCalculatorInstance that would return this instance).  */
	unsigned int numLeaves; /**< The number of leaves in the tree. Set in #createLikelihoodCalcInstance by the argument of the same name */
	unsigned long numPatterns; /**< The number of data patterns. Set in #createLikelihoodCalcInstance by the argument of the same name */
	double * patternWeights; /**< The 0L or a pointer to an array of length `numPatterns` that holds the weight for each pattern. Set in #createLikelihoodCalcInstance by the argument of the same name */
	unsigned int numStates; /**< The number of states in the model. Set in #createLikelihoodCalcInstance by the argument of the same name */
	unsigned int numStateCodeArrays; /**< The number of of data structs to allocate for storing leaf data. Set in #createLikelihoodCalcInstance by the argument of the same name */
	unsigned int numPartialStructs; /**< The f data structures to allocate to store partial likelihood arrays. Set in #createLikelihoodCalcInstance by the argument of the same name */
	unsigned int numProbMats; /**< The . Set in #createLikelihoodCalcInstance by the argument of the same name */
	unsigned int numRescalingsMultipliers; /**< The . Set in #createLikelihoodCalcInstance by the argument of the same name */
	int resourceIndex; /**< The . Set in #createLikelihoodCalcInstance by the argument of the same name */
	long resourcePref; /**< The . Set in #createLikelihoodCalcInstance by the argument of the same name */
	long resourceReq; /**< The . Set in #createLikelihoodCalcInstance by the argument of the same name */
	int beagleInstanceCreated; /**< The . Set in #createLikelihoodCalcInstance by the argument of the same name */


	unsigned int numInstRateModels; /**< The number of "slots" for DSCTModelObj objects. Set in #createLikelihoodCalcInstance by the argument of the same name */
	DSCTModelObj ** probModelArray; /**< The . Set in #createLikelihoodCalcInstance by the argument of the same name */
	const ASRVObj ** asrvAliasForEachModel; /**< The . Length == numInstRateModels */
    
	unsigned int numEigenStorage;
	EigenSolutionStruct ** eigenSolutionStructs;
    
    double * edgeLenScratch; /* length numProbMats*/
    int * probMatIndexScratch; /* length numProbMats */
    double ** probMatScratch; /* numStates x numStates */
    int * stateCodeArrayScratch; /* */
    BeagleOperation * opScratch; /* numPartialStructs long */
    int * waitPartialIndexScratch; /* numPartialStructs long */
    int * categWeightIndexScratch; /* numEigenStorage long */
    double * categWeightScratch; /* numEigenStorage long */
    int * rootPartialIndexScratch; /* numEigenStorage long */
    int * stateFreqIndexScratch; /* numEigenStorage long */
    int * rootRescalerIndexScratch; /* numEigenStorage long */

    
    
};


const struct LikeCalculatorInstance * getLikeCalculatorInstance(long handle);
int calcPrMatsForLCI(const struct LikeCalculatorInstance * lci, int eigenIndex, unsigned numToCalc, const double * edgeLenArray, const int * probMatIndexArray);
int setQMatForLCI(const struct LikeCalculatorInstance * lci, int probMatIndex, const double ** newQMat);

#if defined(API_TRACE_PRINTING) && API_TRACE_PRINTING
        int getTraceModeModelIndex(struct LikeCalculatorInstance * LCI, DSCTModelObj * m);
#endif 


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
