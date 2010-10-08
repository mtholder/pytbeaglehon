#if ! defined(C_PHY_PROB_DEFS_H)
#define C_PHY_PROB_DEFS_H

#ifdef __cplusplus
extern "C"
{
#endif

#if defined(HAVE_INLINE) && HAVE_INLINE
#	define INLINE inline
#else
#	define INLINE
#endif

#if defined(BUILDING_FOR_PYTHON)
#	include "Python.h"
#else
#	define PyObject_HEAD char pyobjectHead;
#	define PyExc_ValueError 1
#	define PyExc_RuntimeError 2
#	define PyExc_IndexError 3
#	define PyExc_TypeError 4
#	define PyObject_New(a,b) ((a*) malloc(sizeof(a)));
#	define PyObject_Del(a) (free((a)));
#	define Py_DECREF(a)
#	define Py_INCREF(a)

	void PyErr_NoMemory();
	void PyErr_SetString(int, const char *);
#endif



#if defined(DEBUG_PRINTING) && DEBUG_PRINTING

#	include <stdio.h>
#	define PYTBEAGLEHON_DEBUG_PRINTF(v) (printf(v))
#	define PYTBEAGLEHON_DEBUG_PRINTF1(v,a) (printf(v,a))
#	define PYTBEAGLEHON_DEBUG_PRINTF2(v,a,aa) (printf(v,a,aa))
#	define PYTBEAGLEHON_DEBUG_PRINTF3(v,a,aa,aaa) (printf(v,a,aa,aaa))
#else

#	define PYTBEAGLEHON_DEBUG_PRINTF(v)
#	define PYTBEAGLEHON_DEBUG_PRINTF1(v,a)
#	define PYTBEAGLEHON_DEBUG_PRINTF2(v,a,aa)
#	define PYTBEAGLEHON_DEBUG_PRINTF3(v,a,aa,aaa)

#endif

#include <assert.h>
#define PYTBEAGLE_ASSERT(x) assert(x);




/*


struct LikeCalculatorInstance;

int likeCalcInstanceValid(struct LikeCalculatorInstance * b);
int getProbMatIndex(struct LikeCalculatorInstance * b, int *indices, int numIndices);
*/






#ifdef __cplusplus
}
/* extern "C" */
#endif

#endif
