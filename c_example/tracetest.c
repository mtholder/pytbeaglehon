#include "pytbeaglehon/ccore/asrv.h"
#include "pytbeaglehon/ccore/calc_instance.h"
#include "pytbeaglehon/ccore/discrete_state_model.h"
#include "pytbeaglehon/ccore/internal_like_calc_env.h"
#include <stdlib.h>
#include <stdio.h>

int eigenIndex, numToCalc;
DSCTModelObj * dsct_model_obj = 0L;
double lnL;
int main(int argc, char * argv[]) {
/* cAPI Call */ int numLeaves=4; long numPatterns=16 ; const double * patternWeights=0L; int numStates=4;int numStateCodeArrays=4; int numPartialStructs=6; int numInstRateModels=1 ; const ASRVObj ** asrvObjectArray=0L;int numProbMats=14; int numEigenStorage=2; int numRescalingsMultipliers=4; int resourceIndex=0; long resourcePref=524736; long resourceReq=524736;
/* cAPI Call */ long handle = createLikelihoodCalcInstance(numLeaves, numPatterns, patternWeights, numStates, numStateCodeArrays, numPartialStructs, numInstRateModels, asrvObjectArray, numProbMats, numEigenStorage, numRescalingsMultipliers, resourceIndex, resourcePref, resourceReq);
/* cAPI Call */ struct LikeCalculatorInstance * LCI = getLikeCalculatorInstance(0);
/* cAPI Call */ handle = 0; int scarr0[] = {0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3};
/* cAPI Call */ setStateCodeArray(handle, 0, scarr0);
/* cAPI Call */ handle = 0; int scarr1[] = {0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3};
/* cAPI Call */ setStateCodeArray(handle, 1, scarr1);
/* cAPI Call */ handle = 0; int scarr2[] = {0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3};
/* cAPI Call */ setStateCodeArray(handle, 2, scarr2);
/* cAPI Call */ handle = 0; int scarr3[] = {0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3};
/* cAPI Call */ setStateCodeArray(handle, 3, scarr3);
/* cAPI Call */ dsct_model_obj = LCI->probModelArray[0];
/* cAPI Call */ dsct_model_obj->qMat[0][0] = -0.862245; dsct_model_obj->qMat[0][1] = 0.125000; dsct_model_obj->qMat[0][2] = 0.612245; dsct_model_obj->qMat[0][3] = 0.125000;
/* cAPI Call */ dsct_model_obj->eigenCalcIsDirty = 1;
/* cAPI Call */ dsct_model_obj = LCI->probModelArray[0];
/* cAPI Call */ dsct_model_obj->eigenBufferIndex=0; dsct_model_obj->eigenCalcIsDirty = 1; recalc_eigen_mat(dsct_model_obj);
/* cAPI Call */ handle=0; eigenIndex=0; numToCalc=1;LCI->edgeLenScratch[0] = 0.000000; ;
/* cAPI Call */ LCI->probMatIndexScratch[0] = 0; 
/* cAPI Call */ calcPrMats(handle, eigenIndex, numToCalc, LCI->edgeLenScratch, LCI->probMatIndexScratch);
/* cAPI Call */ handle=0; eigenIndex=0; numToCalc=1;LCI->edgeLenScratch[0] = 0.000000; ;
/* cAPI Call */ LCI->probMatIndexScratch[0] = 8; 
/* cAPI Call */ calcPrMats(handle, eigenIndex, numToCalc, LCI->edgeLenScratch, LCI->probMatIndexScratch);
/* cAPI Call */ handle=0; eigenIndex=0; numToCalc=1;LCI->edgeLenScratch[0] = 0.000000; ;
/* cAPI Call */ LCI->probMatIndexScratch[0] = 6; 
/* cAPI Call */ calcPrMats(handle, eigenIndex, numToCalc, LCI->edgeLenScratch, LCI->probMatIndexScratch);
/* cAPI Call */ handle=0; eigenIndex=0; numToCalc=1;LCI->edgeLenScratch[0] = 0.000000; ;
/* cAPI Call */ LCI->probMatIndexScratch[0] = 11; 
/* cAPI Call */ calcPrMats(handle, eigenIndex, numToCalc, LCI->edgeLenScratch, LCI->probMatIndexScratch);
/* cAPI Call */ handle=0; eigenIndex=0; numToCalc=1;LCI->edgeLenScratch[0] = 0.000000; ;
/* cAPI Call */ LCI->probMatIndexScratch[0] = 4; 
/* cAPI Call */ calcPrMats(handle, eigenIndex, numToCalc, LCI->edgeLenScratch, LCI->probMatIndexScratch);
/* cAPI Call */ handle=0; eigenIndex=0; numToCalc=1;LCI->edgeLenScratch[0] = 0.010000; ;
/* cAPI Call */ LCI->probMatIndexScratch[0] = 2; 
/* cAPI Call */ calcPrMats(handle, eigenIndex, numToCalc, LCI->edgeLenScratch, LCI->probMatIndexScratch);
/* cAPI Call */ handle = 0; LCI->opScratch[0] = partialOperation(6, -1, -1, 2, 0, 3, 8); LCI->opScratch[1] = partialOperation(4, -1, -1, 0, 6, 1, 11); LCI->opScratch[2] = partialOperation(8, -1, -1, 4, 4, 6, 2); 
/* cAPI Call */ LCI->waitPartialIndexScratch[0] = 8; 
/* cAPI Call */ calcPartials(handle, LCI->opScratch, 3, LCI->waitPartialIndexScratch, 1);
/* cAPI Call */ handle = 0; LCI->categWeightIndexScratch[0] = 0; LCI->categWeightScratch[0] = 1.000000; 
/* cAPI Call */ setSingletonCategoryWeights(handle, LCI->categWeightIndexScratch, LCI->categWeightScratch, 1);
/* cAPI Call */ handle = 0; eigenIndex = 0; LCI->categWeightScratch[0] = 0.300000; LCI->categWeightScratch[1] = 0.250000; LCI->categWeightScratch[2] = 0.200000; LCI->categWeightScratch[3] = 0.250000; 
/* cAPI Call */ setStateFreq(handle, eigenIndex, LCI->categWeightScratch);
/* cAPI Call */ dsct_model_obj = LCI->probModelArray[0];
/* cAPI Call */ dsct_model_obj->qMat[0][0] = -0.862245; dsct_model_obj->qMat[0][1] = 0.125000; dsct_model_obj->qMat[0][2] = 0.612245; dsct_model_obj->qMat[0][3] = 0.125000;
/* cAPI Call */ dsct_model_obj->eigenCalcIsDirty = 1;
/* cAPI Call */ dsct_model_obj = LCI->probModelArray[0];
/* cAPI Call */ dsct_model_obj->eigenBufferIndex=1; dsct_model_obj->eigenCalcIsDirty = 1; recalc_eigen_mat(dsct_model_obj);
/* cAPI Call */ LCI->rootPartialIndexScratch[0] = 8; /* cAPI Call */ LCI->categWeightIndexScratch[0] = 1; /* cAPI Call */ LCI->stateFreqIndexScratch[0] = 1; /* cAPI Call */ LCI->rootRescalerIndexScratch[0] = 4; 
/* cAPI Call */ calcRootLnL(handle, LCI->rootPartialIndexScratch, LCI->categWeightIndexScratch, LCI->stateFreqIndexScratch, LCI->rootRescalerIndexScratch, 1, &lnL);
return 0;}
