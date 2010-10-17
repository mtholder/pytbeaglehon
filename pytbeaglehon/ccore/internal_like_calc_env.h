/**  Provides a handle for the memory allocated for likelihood calculations (other
    than the ASRV structurs
 */
#if ! defined(INTERNAL_LIKE_CALC_ENV_H)
#define INTERNAL_LIKE_CALC_ENV_H
#ifdef __cplusplus
extern "C" 
{
#endif

#include "discrete_state_model.h"


/* Each LikeCalculatorInstance contains the internal storage needed to
    calculate likelihoods for one subset in a (possibly) partitioned matrix.

    Multiple Models may be associate with the instance, but all model share the
    same number of states.

*/
struct LikeCalculatorInstance {
    int beagleInstanceIndex;

	unsigned int numLeaves;
	unsigned long numPatterns;
	double * patternWeights;
	unsigned int numStates;
	unsigned int numRateCategories;
	unsigned int numStateCodeArrays;
	unsigned int numPartialStructs;
	unsigned int numProbMats;
	unsigned int numRescalingsMultipliers;
	int resourceIndex;
	long resourcePref;
	long resourceReq;
	int beagleInstanceCreated;


	unsigned int numInstRateModels;
	DSCTModelObj ** probModelArray;
	const ASRVObj ** asrvAliasForEachModel; /* length = numInstRateModels */
    
	unsigned int numEigenStorage;
	EigenSolutionStruct ** eigenSolutionStructs;
    
    double * edgeLenScratch; /* length numProbMats*/
    int * probMatIndexScratch; /* length numProbMats */
    
};


struct LikeCalculatorInstance * getLikeCalculatorInstance(long handle);

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
