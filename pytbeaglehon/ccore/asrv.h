/**  Among-site rate variation calculation in C.  Does not use BEAGLE.
 */
#if ! defined(ASRV_H)
#define ASRV_H
#ifdef __cplusplus
extern "C" 
{
#endif




#include "pytbeaglehon_defs.h"


/**
 * stores an array of doubles and an array of frequencies representing the
 * a probability for each element in the array.
 * An example usage is for storing an array of rates for n_categories and the
 * probability that a site would belong to each of the rate categories.
 */
typedef struct {
	PyObject_HEAD
	unsigned n;
	int style; /* enum facet for to indicate type -- this is a hack
				  0 ="use median for rate categ"
				  1 ="use mean for rate categ"
			   */
	double param; /*shape param of the gamma distribution -- this is a hack*/
	double * val;
	double * freq;
} ASRVObj;



ASRVObj* asrv_obj_new(unsigned dim, int style, double param);
void asrv_obj_dtor(ASRVObj* asrh);
void internal_asrv_set_shape(ASRVObj *asrh, double val);






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
