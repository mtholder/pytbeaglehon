/**
  adaptor code for filling in ASRVObj

    This file does not actually use python-c extension (calls like PyObject_Del are
        defined to free in pytbeaglehon_defs.h if BUILDING_FOR_PYTHON is not defined).
*/

#include "asrv.h"
#include "phylo_util.h"

void internal_asrv_set_shape(ASRVObj *asrh, double val) {
	PYTBEAGLEHON_DEBUG_PRINTF("In internal_asrv_set_shape\n");
	double beta;
	asrh->param = val;
	beta = 1.0/val;
	DiscreteGamma(asrh->freq, asrh->val, val, beta, asrh->n, asrh->style);
}



void asrv_obj_dtor(ASRVObj * asrh) {
	PYTBEAGLEHON_DEBUG_PRINTF("In asrh dtor\n");
	if (asrh == 0L)
		return;
	if (asrh->val != 0L)
		free(asrh->val);
	if (asrh->freq != 0L)
		free(asrh->freq);
	PyObject_Del(asrh);
}





/**
  pytbeaglehon phylogenetic likelihood caluclations using beaglelib.

  Copyright 2010 Mark T. Holder (except where noted). All rights reserved.

  See "LICENSE.txt" for terms and conditions of usage.
*/
