/**
    Storage structures for a model of evolution (and eigensolutions).
 */
#if ! defined(DSCT_MODEL_H)
#define DSCT_MODEL_H

#ifdef __cplusplus
extern "C" 
{
#endif

#include "pytbeaglehon_defs.h"

#include "calc_instance.h"

typedef struct {
		double ** workMat; /*dim by dim*/
		double * dWork;  /*len dim*/
		int * iWork; /*len dim*/	
} EigenCalcScratchpad;

typedef struct {
    int beagleEigenBufferIndex;
	double * eigenValues; /* alias */
	double * imEigenValues; /* alias */
	double ** eigenVectors;
	double ** invEigenVectors;
	EigenCalcScratchpad * scratchPad;
	unsigned int dim;
	double *** matDelHandle; 
	double * arrDelHandle; 
} EigenSolutionStruct;

EigenSolutionStruct * eigenSolutionStructNew(unsigned dim);
void eigenSolutionStructDtor(EigenSolutionStruct * p);

struct LikeCalculatorInstance;
/* Define the type that corresponds to a model (holds the Q-Matrix and other
	temporary fields used to speed up calculation).
*/
typedef struct {
	PyObject_HEAD
	unsigned dim;
	double **qMat;
	struct LikeCalculatorInstance * likeCalcInstanceAlias;
	int eigenBufferIndex; /* index where this eigen system is currently stored */
	int eigenCalcIsDirty; /*0 if the eigen calculations are up-to-date*/
} DSCTModelObj;

DSCTModelObj * dsctModelNew(unsigned dim);

const DSCTModelObj ** getModelList(long instanceHandle, unsigned int *numModels);

int set_eigen_index_for_model(DSCTModelObj * dsct_model_obj, int eigenIndex);
int calc_eigen_mat(DSCTModelObj * dsct_model_obj, int eigenIndex);
int recalc_eigen_mat(DSCTModelObj *mod);
int set_q_mat_for_model(DSCTModelObj * dsct_model_obj, const double ** newQMat);



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
