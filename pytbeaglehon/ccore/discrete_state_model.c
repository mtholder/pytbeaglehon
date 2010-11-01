#include "internal_like_calc_env.h"
#include "phylo_util.h"
#include <assert.h>
#include <stdlib.h>
#include <libhmsbeagle/beagle.h>




EigenCalcScratchpad * eigenSolutionScratchpadNew(unsigned dim);
void eigenSolutionScratchpadDtor(EigenCalcScratchpad *);
void cdsctm_dtor(DSCTModelObj* dsct_model);


#if defined(BUILDING_FOR_PYTHON)

	PyTypeObject dsct_model_type;

	PyTypeObject dsct_model_type = {
	PyObject_HEAD_INIT(0)	  /* initialize to 0 to ensure Win32 portability  */
	0,						  /* ob_size */
	"dsct_model",					/* tp_name */
	sizeof(DSCTModelObj),		/* tp_basicsize */
	0,						  /* tp_itemsize */
	/* methods */
	(destructor)cdsctm_dtor, /* tp_dealloc */
	/* implied by ISO C: all zeros thereafter, i.e., no other method */
	};


#endif



DSCTModelObj * dsctModelNew(unsigned dim)  {
	assert(dim > 1);
	DSCTModelObj * dsct_model = PyObject_New(DSCTModelObj, &dsct_model_type);
	PYTBEAGLEHON_DEBUG_PRINTF1("In dsctModelNew. Allocated address = %ld\n", (long) dsct_model);


	if (dsct_model) {
		dsct_model->dim = dim;
		dsct_model->eigenBufferIndex = -1;
		dsct_model->eigenCalcIsDirty = 1;
		dsct_model->likeCalcInstanceAlias = 0L;
		dsct_model->qMat = 0L;
		dsct_model->qMat = allocateDblMatrix(dim, dim);
		if (dsct_model->qMat == 0L)
			goto errorExit;
	}
	return dsct_model;
	errorExit:
		PYTBEAGLEHON_DEBUG_PRINTF("In cdsctm_ctor errorExit\n");
		cdsctm_dtor(dsct_model);
		return 0L;
}


void cdsctm_dtor(DSCTModelObj* dsct_model) {
	PYTBEAGLEHON_DEBUG_PRINTF1("In cdsctm_dtor address = %ld\n", (long) dsct_model);
	if (dsct_model == 0L)
		return;
	if (dsct_model->qMat)
		freeDblMatrix(dsct_model->qMat);
	PyObject_Del(dsct_model);
}



EigenSolutionStruct * eigenSolutionStructNew(unsigned dim) {
    assert(dim > 1);
	EigenSolutionStruct * essPtr = (EigenSolutionStruct *)malloc(sizeof(EigenSolutionStruct));
	if (essPtr == 0L) {
		PYTBEAGLEHON_DEBUG_PRINTF("could not alloc EigenSolutionStruct\n");
		return 0L;
	}

	essPtr->scratchPad = 0L;
	essPtr->dim = dim;
	essPtr->eigenVectors = 0L;
	essPtr->invEigenVectors = 0L;
	essPtr->eigenValues = 0L;
	essPtr->imEigenValues = 0L;

    essPtr->beagleEigenBufferIndex = -1;
	unsigned int lenArrDelHandle = 2*dim;
	essPtr->matDelHandle = 0L;
	essPtr->arrDelHandle = 0L;

	essPtr->matDelHandle = allocateDbl3DMatrix(2, dim, dim);
	if (essPtr->matDelHandle == 0L) {
		PYTBEAGLEHON_DEBUG_PRINTF("alloc fail in EigenSolutionStruct->matDelHandle\n");
		goto errorExit;
	}
	essPtr->arrDelHandle = (double*)malloc(lenArrDelHandle*sizeof(double));
	if (essPtr->arrDelHandle == 0L) {
		PYTBEAGLEHON_DEBUG_PRINTF("alloc fail in EigenSolutionStruct->arrDelHandle\n");
		goto errorExit;
	}
	essPtr->scratchPad = eigenSolutionScratchpadNew(dim);
	if (essPtr->scratchPad == 0L){
		PYTBEAGLEHON_DEBUG_PRINTF("alloc fail in EigenSolutionStruct->scratchPad\n");
		goto errorExit;
	}
	essPtr->eigenVectors = essPtr->matDelHandle[0];
	essPtr->invEigenVectors = essPtr->matDelHandle[1];
	essPtr->eigenValues = essPtr->arrDelHandle;
	essPtr->imEigenValues = essPtr->eigenValues + dim;
    PYTBEAGLEHON_DEBUG_PRINTF1("Allocated eigenSolutionStructNew at %ld\n", (long) essPtr);
	return essPtr;

	errorExit:
		PYTBEAGLEHON_DEBUG_PRINTF("eigenSolutionStructNew errorExit\n");
		eigenSolutionStructDtor(essPtr);
		return 0L;

}

/* Used internally in the lib to allocate the workspace for eigensolution calculations */
EigenCalcScratchpad * eigenSolutionScratchpadNew(unsigned dim) {
	EigenCalcScratchpad * essPtr = (EigenCalcScratchpad *)malloc(sizeof(EigenCalcScratchpad));
	if (essPtr == 0L) {
		PYTBEAGLEHON_DEBUG_PRINTF("EigenCalcScratchpad errorExit\n");
		return 0L;
	}
	essPtr->workMat = 0L;
	essPtr->dWork = 0L;
	essPtr->iWork = 0L;
	essPtr->iWork = (int *)malloc(dim*sizeof(int));
	if (essPtr->iWork == 0L) {
		PYTBEAGLEHON_DEBUG_PRINTF("alloc fail in EigenCalcScratchpad->iWork\n");
		goto errorExit;
	}
	essPtr->dWork = (double *)malloc(dim*sizeof(double));
	if (essPtr->dWork == 0L) {
		PYTBEAGLEHON_DEBUG_PRINTF("alloc fail in EigenCalcScratchpad->dWork\n");
		goto errorExit;
	}
	essPtr->workMat = allocateDblMatrix(dim, dim);
	if (essPtr->workMat == 0L) {
		PYTBEAGLEHON_DEBUG_PRINTF("alloc fail in EigenCalcScratchpad->workMat\n");
		goto errorExit;
	}
    PYTBEAGLEHON_DEBUG_PRINTF4("Allocated EigenCalcScratchpad at %ld ( iwork=%ld, dwork=%ld, workMat=%ld )\n", (long) essPtr, 
                                                                                                               (long) essPtr->iWork,
                                                                                                               (long) essPtr->dWork,
                                                                                                               (long) essPtr->workMat);
	return essPtr;

	errorExit:
		PYTBEAGLEHON_DEBUG_PRINTF("EigenCalcScratchpad errorExit\n");
		eigenSolutionScratchpadDtor(essPtr);
		return 0L;
}



void eigenSolutionScratchpadDtor(EigenCalcScratchpad * p) {
	if (p == 0L)
		return;
    PYTBEAGLEHON_DEBUG_PRINTF4("Deleting EigenCalcScratchpad at %ld ( iwork=%ld, dwork=%ld, workMat=%ld )\n", (long) p, 
                                                                                                               (long) p->iWork,
                                                                                                               (long) p->dWork,
                                                                                                               (long) p->workMat);
	if (p->workMat != 0L) {
		freeDblMatrix(p->workMat);
	}
	if (p->dWork != 0L) {
		free(p->dWork);
	}
	if (p->iWork != 0L) {
		free(p->iWork);
	}
	free(p);
    PYTBEAGLEHON_DEBUG_PRINTF("finished scratchpad dtor\n"); 
}
void eigenSolutionStructDtor(EigenSolutionStruct * p) {
	if (p == 0L)
		return;
    PYTBEAGLEHON_DEBUG_PRINTF1("Deleting EigenSolutionStruct at %ld\n", (long) p);
	if (p->matDelHandle)
		freeDbl3DMatrix(p->matDelHandle);
	if (p->arrDelHandle)
		free(p->arrDelHandle);
	eigenSolutionScratchpadDtor(p->scratchPad);
	free(p);
    PYTBEAGLEHON_DEBUG_PRINTF("finished eigenSolutionStructDtor\n"); 
}



EigenSolutionStruct * getEigenSolutionStruct(DSCTModelObj * mod) {
	assert(mod);
	struct LikeCalculatorInstance * lce = mod->likeCalcInstanceAlias;
	if (!mod && !(lce))
	    return 0L;
	
	if (lce->numEigenStorage <= mod->eigenBufferIndex)
	    return 0L;
	return lce->eigenSolutionStructs[mod->eigenBufferIndex];
}



/**
 * Recalculates the eigensystem and temporaries.  Returns 0 if any steps fail
 *  or the complex eigenvalues are encountered.
 *
 *	\Returns 0 to indicate failure (if this happens PyErr_SetString will have
 *		been called).
 */
int recalc_eigen_mat(DSCTModelObj *mod) {
	unsigned is_complex;
	int rc;
	struct LikeCalculatorInstance * lce = mod->likeCalcInstanceAlias;
	EigenSolutionStruct * eigenSolutionStruct = getEigenSolutionStruct(mod);
	EigenCalcScratchpad * eigenCalcScratchpad = (eigenSolutionStruct == 0L ? 0L : eigenSolutionStruct->scratchPad);

    if ((eigenSolutionStruct == 0L) || (eigenCalcScratchpad == 0L)) {
		PyErr_SetString(PyExc_RuntimeError, "Could not get reference to eigen solution storage .");
		return 0;
	}
#	if defined(PRINTING_LOTS) && PRINTING_LOTS
		printQMat(mod->qMat, mod->dim);
#	endif

    PYTBEAGLEHON_DEBUG_PRINTF2("Calculating eigensolution for model at %ld in eigenstorage %d\n", (long)mod, eigenSolutionStruct->beagleEigenBufferIndex);

	rc = get_eigens(mod->dim,
					mod->qMat,
					eigenSolutionStruct->eigenValues,
					eigenSolutionStruct->imEigenValues,
					eigenSolutionStruct->eigenVectors,
					eigenSolutionStruct->invEigenVectors,
					eigenCalcScratchpad->workMat,
					eigenCalcScratchpad->dWork,
					eigenCalcScratchpad->iWork,
					&is_complex);
	if (rc == 0) {
		if (is_complex == 1)
			PyErr_SetString(PyExc_RuntimeError, "Error: Complex eigenvalues found.");
		return 0;
	}
	
#   if defined(BEAGLE_API_TRACE_PRINTING) && BEAGLE_API_TRACE_PRINTING
        int rowInd, colInd;
        for (rowInd = 0; rowInd < lce->numStates; ++rowInd) {
            PYTBEAGLEHON_DEBUG_PRINTF3("/* BEAGLE_API Call */ inst%dEigenVal[%d] = %lf;\n", lce->beagleInstanceIndex, rowInd, eigenSolutionStruct->eigenValues[rowInd]);
            for (colInd = 0; colInd < lce->numStates; ++colInd) {
                PYTBEAGLEHON_DEBUG_PRINTF4("/* BEAGLE_API Call */ inst%dEigenVec[%d][%d] = %lf;\n", lce->beagleInstanceIndex, rowInd, colInd, eigenSolutionStruct->eigenVectors[rowInd][colInd]);
                PYTBEAGLEHON_DEBUG_PRINTF4("/* BEAGLE_API Call */ inst%dInvEigenVec[%d][%d] = %lf;\n", lce->beagleInstanceIndex, rowInd, colInd, eigenSolutionStruct->invEigenVectors[rowInd][colInd]);
            }
        }
        PYTBEAGLEHON_DEBUG_PRINTF5("/* BEAGLE_API Call */ beagleSetEigenDecomposition(%d, %d, inst%dEigenVec[0], inst%dInvEigenVec[0], inst%dEigenVal);\n", lce->beagleInstanceIndex, eigenSolutionStruct->beagleEigenBufferIndex, lce->beagleInstanceIndex, lce->beagleInstanceIndex, lce->beagleInstanceIndex);
#   endif
	rc = beagleSetEigenDecomposition(lce->beagleInstanceIndex,
                                     eigenSolutionStruct->beagleEigenBufferIndex,
                                     eigenSolutionStruct->eigenVectors[0],
                                     eigenSolutionStruct->invEigenVectors[0],
                                     eigenSolutionStruct->eigenValues);
    if (rc != BEAGLE_SUCCESS) {
		PyErr_SetString(PyExc_RuntimeError, "Error: call to beagleSetEigenDecomposition failed.");
		return 0;
	}
	mod->eigenCalcIsDirty = 0;
	return 1;
}



/**
  pytbeaglehon phylogenetic likelihood caluclations using beaglelib.

  Copyright 2010 Mark T. Holder (except where noted). All rights reserved.

  See "LICENSE.txt" for terms and conditions of usage.
*/
 
