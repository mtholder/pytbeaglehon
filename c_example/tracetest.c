#include "pytbeaglehon/ccore/asrv.h"
#include "pytbeaglehon/ccore/calc_instance.h"
#include "pytbeaglehon/ccore/discrete_state_model.h"
#include <stdlib.h>
#include <stdio.h>

int eigenIndex, numToCalc;
double * edgeLen;
double * categWeightScratch;
int * categWeightIndexScratch;
int * probMatIndexScratch;
int * waitPartialIndexScratch;
BeagleOperation * opScratch;

int main(int argc, char * argv[]) {
    edgeLen = (double *)malloc(1000*sizeof(double));
    categWeightScratch = (double *)malloc(1000*sizeof(double));
    categWeightIndexScratch = (int *)malloc(1000*sizeof(int));
    probMatIndexScratch = (int *)malloc(1000*sizeof(int));
    waitPartialIndexScratch = (int *)malloc(1000*sizeof(int));
    opScratch = (BeagleOperation *)malloc(1000*sizeof(BeagleOperation));
/* cAPI Call */ int numLeaves=4; long numPatterns=16 ; const double * patternWeights=0L; int numStates=4;int numStateCodeArrays=4; int numPartialStructs=6; int numInstRateModels=1 ; const ASRVObj ** asrvObjectArray=0L;int numProbMats=14; int numEigenStorage=2; int numRescalingsMultipliers=4; int resourceIndex=0; long resourcePref=524736; long resourceReq=524736;
/* cAPI Call */ long handle = createLikelihoodCalcInstance(numLeaves, numPatterns, patternWeights, numStates, numStateCodeArrays, numPartialStructs, numInstRateModels, asrvObjectArray, numProbMats, numEigenStorage, numRescalingsMultipliers, resourceIndex, resourcePref, resourceReq);
/* cAPI Call */ DSCTModelObj * dsct_model_obj = 0L;
/* cAPI Call */ DSCTModelObj * modObj0 = dsctModelNew(4);
/* cAPI Call */ handle = 0; int scarr0[] = {0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3};
/* cAPI Call */ setStateCodeArray(handle, 0, scarr0);
/* cAPI Call */ handle = 0; int scarr1[] = {0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3};
/* cAPI Call */ setStateCodeArray(handle, 1, scarr1);
/* cAPI Call */ handle = 0; int scarr2[] = {0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3};
/* cAPI Call */ setStateCodeArray(handle, 2, scarr2);
/* cAPI Call */ handle = 0; int scarr3[] = {0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3};
/* cAPI Call */ setStateCodeArray(handle, 3, scarr3);
/* cAPI Call */ dsct_model_obj = modObj0;
/* cAPI Call */ dsct_model_obj->qMat[0][0] = -0.862245; dsct_model_obj->qMat[0][1] = 0.125000; dsct_model_obj->qMat[0][2] = 0.612245; dsct_model_obj->qMat[0][3] = 0.125000;
/* cAPI Call */ dsct_model_obj->eigenCalcIsDirty = 1;
/* cAPI Call */ dsct_model_obj->eigenBufferIndex=0; dsct_model_obj->eigenCalcIsDirty = 1; recalc_eigen_mat(dsct_model_obj);
/* cAPI Call */ handle=0; eigenIndex=0; numToCalc=1;edgeLen[0] = 0.000000; ;
/* cAPI Call */ probMatIndexScratch[0] = 6; 
/* cAPI Call */ calcPrMats(handle, eigenIndex, numToCalc, edgeLen, probMatIndexScratch);
/* cAPI Call */ handle=0; eigenIndex=0; numToCalc=1;edgeLen[0] = 0.000000; ;
/* cAPI Call */ probMatIndexScratch[0] = 8; 
/* cAPI Call */ calcPrMats(handle, eigenIndex, numToCalc, edgeLen, probMatIndexScratch);
/* cAPI Call */ handle=0; eigenIndex=0; numToCalc=1;edgeLen[0] = 0.000000; ;
/* cAPI Call */ probMatIndexScratch[0] = 10; 
/* cAPI Call */ calcPrMats(handle, eigenIndex, numToCalc, edgeLen, probMatIndexScratch);
/* cAPI Call */ handle=0; eigenIndex=0; numToCalc=1;edgeLen[0] = 0.000000; ;
/* cAPI Call */ probMatIndexScratch[0] = 11; 
/* cAPI Call */ calcPrMats(handle, eigenIndex, numToCalc, edgeLen, probMatIndexScratch);
/* cAPI Call */ handle=0; eigenIndex=0; numToCalc=1;edgeLen[0] = 0.000000; ;
/* cAPI Call */ probMatIndexScratch[0] = 4; 
/* cAPI Call */ calcPrMats(handle, eigenIndex, numToCalc, edgeLen, probMatIndexScratch);
/* cAPI Call */ handle=0; eigenIndex=0; numToCalc=1;edgeLen[0] = 0.010000; ;
/* cAPI Call */ probMatIndexScratch[0] = 2; 
/* cAPI Call */ calcPrMats(handle, eigenIndex, numToCalc, edgeLen, probMatIndexScratch);
/* cAPI Call */ handle = 0; opScratch[0] = partialOperation(6, -1, -1, 2, 6, 3, 8); opScratch[1] = partialOperation(8, -1, -1, 0, 10, 1, 11); opScratch[2] = partialOperation(9, -1, -1, 8, 4, 6, 2); 
/* cAPI Call */ waitPartialIndexScratch[0] = 9; 
/* cAPI Call */ calcPartials(handle, opScratch, 3, waitPartialIndexScratch, 1);
/* cAPI Call */ handle = 0; categWeightIndexScratch[0] = 0; categWeightScratch[0] = 1.000000; 
/* cAPI Call */ setSingletonCategoryWeights(handle, categWeightIndexScratch, categWeightScratch, 1);
/* cAPI Call */ dsct_model_obj = modObj0;
/* cAPI Call */ dsct_model_obj->qMat[0][0] = -0.862245; dsct_model_obj->qMat[0][1] = 0.125000; dsct_model_obj->qMat[0][2] = 0.612245; dsct_model_obj->qMat[0][3] = 0.125000;
/* cAPI Call */ dsct_model_obj->eigenCalcIsDirty = 1;
/* cAPI Call */ dsct_model_obj->eigenBufferIndex=1; dsct_model_obj->eigenCalcIsDirty = 1; recalc_eigen_mat(dsct_model_obj);
 return 0; }
