#include "pytbeaglehon_defs.h"
#include "phylo_util.h"
#include <float.h>
#include <math.h>
#include <stdlib.h>

/**
  pytbeaglehon phylogenetic likelihood caluclations using beaglelib.

  Copyright 2010 Mark T. Holder (except where noted). 
  All rights reserved.
  The DiscreteGamma function declared here is implemented using code from
    Ziheng Yang's PAML package

  The get_eigens function declared here is implemented using code from
    MrBayes (Huelsenbeck, Ronquist, van der Mark)

  See "LICENSE.txt" for terms and conditions of usage.
  

*/

void ** mallocZeroedPointerArray(unsigned i) {
	unsigned j;
	void ** p = (void**) malloc(i*sizeof(void *));
	if (p) {
		for (j = 0; j < i; ++j)
			p[j] = 0L;
	}
	return p;
}



static void CopyDoubleMatrices (int dim, const double **from, double **to);
static int compute_eigen_system (int dim, double **a, double *v, double *vi, double **u, int *iwork, double *dwork, unsigned *is_complex);
static int InvertMatrix (int dim,
						 double **a, /*dim by dim - valid on input, but overwritten */
						 double *dWork, /* scratch array, len dim */
						 int *iWork, /* scratch array, len dim */
						 double **aInv/*dim by dim - valid on completion unless 0 is returned */
						 );
static void LUBackSubstitution (int dim, double **a, int *iWork, double *b);
static int LUDecompose (int dim, double **a, double *dWork, int *iWork);
static const double MB_ETA_TOL = 1E-30;
static const double MB_TINY_TOL = 1.0e-20;



/**
 * (internal) allocates a 2-D array where the memory for the doubles is a
 * contiguous block.
 *
 *	\assert (n_rows > 0) and (n_cols > 0)
 */
double **allocateDblMatrix(unsigned n_rows, unsigned n_cols) {
	double **mat_p;
	double *arr;
	const int arr_len = n_rows*n_cols*(sizeof(double));
	unsigned i, j;
	PYTBEAGLEHON_DEBUG_PRINTF1("In allocateDblMatrix len = %d \n", n_rows*n_cols);
	assert(n_rows > 0);
	assert(n_cols > 0);
	mat_p = (double**)malloc(n_rows*sizeof(double*));
	if (mat_p == 0L) {
		PyErr_NoMemory();
		return 0L;
	}
	arr = (double*)malloc(arr_len);
	if (arr == 0L) {
		free(mat_p);
		PyErr_NoMemory();
		return 0L;
	}
	for (i = 0; i < n_rows; ++i) {
		mat_p[i] = arr;
		for (j = 0; j < n_cols; ++j)
			*arr++ = 0.0;
	}
	return mat_p;
}
/**
 * (internal) allocates a 2-D array where the memory for the doubles is a
 * contiguous block.
 *
 *	\assert (nMatrices > 0) and (n_rows > 0) and (n_cols > 0)
 */
double ***allocateDbl3DMatrix(unsigned nm, unsigned nr, unsigned nc) {
	double ***mat_list_p;
	double **mat_walk;
	double *arr;
	const int arr_len = nm*nr*nc*(sizeof(double));
	unsigned i, j, k;
	PYTBEAGLEHON_DEBUG_PRINTF1("In allocateDbl3DMatrix len =%d\n", arr_len);
	assert(nm > 0);
	assert(nr > 0);
	assert(nc > 0);
	mat_list_p = (double***)malloc(nm*sizeof(double**));
	if (mat_list_p == 0L) {
		PyErr_NoMemory();
		return 0L;
	}
	mat_walk = (double**)malloc(nr*nm*sizeof(double*));
	if (mat_walk == 0L) {
		free(mat_list_p);
		PyErr_NoMemory();
		return 0L;
	}
	arr = (double*)malloc(arr_len);
	if (arr == 0L) {
		free(mat_list_p);
		free(mat_walk);
		PyErr_NoMemory();
		return 0L;
	}
	for (i = 0; i < nm; ++i) {
		mat_list_p[i] = mat_walk;
		for (j = 0; j < nr; ++j) {
			*mat_walk++ = arr;
			for (k = 0; k < nc; ++k) {
				*arr++ = 0.0;
			}
		}
	}
	return mat_list_p;
}
/**
 * (internal) frees a 2-D array allocated with allocateDblMatrix
 *
 *	\assert (n_rows > 0) and (n_cols > 0)
 */
void freeDblMatrix(double **p) {
	PYTBEAGLEHON_DEBUG_PRINTF("In freeDblMatrix\n");
	if (p == 0L)
		return;
	if (p[0] != 0L)
		free(p[0]);
	free(p);
}
void freeDbl3DMatrix(double ***p) {
	PYTBEAGLEHON_DEBUG_PRINTF("In freeDbl3DMatrix\n");
	if (p == 0L)
		return;
	if (p[0] != 0L) {
		if (p[0][0] != 0L) {
			free(p[0][0]);
		}
		free(p[0]);
	}
	free(p);
}
/* Phylogenetic and Numerical code */
/* from Ziheng Yang's PAML */
/* returns the incomplete gamma ratio I(x,alpha) where x is the upper
		   limit of the integration and alpha is the shape parameter.
   returns (-1) if in error
   ln_gamma_alpha = ln(Gamma(alpha)), is almost redundant.
   (1) series expansion		if (alpha < x || x <= 1)
   (2) continued fraction	otherwise
   RATNEST FORTRAN by
   Bhattacharjee GP (1970) The incomplete gamma integral.  Applied Statistics,
   19: 285-287 (AS32)
*/
double IncompleteGamma (double x, double alpha, double ln_gamma_alpha) {
	int i;
	double p = alpha, g = ln_gamma_alpha;
	const double ACCURACY_TOL = 1e-10;
	const double OVERFLOW_VAL = 1e60;
	double factor, gin = 0, rn = 0, a = 0,b = 0,an = 0,dif = 0, term = 0, pn[6];
	if (x == 0)
		return 0;
	if (x < 0 || p <= 0)
		return -1;
	factor = exp(p*log(x) - x - g);
	if (x > 1 && x >= p) {
		/* continued fraction */
		a = 1 - p;
		b = a + x + 1;
		term = 0;
		pn[0] = 1;
		pn[1] = x;
		pn[2] = x + 1;
		pn[3] = x*b;
		gin = pn[2]/pn[3];
		for (;;){
			a++;
			b += 2;
			term++;
			an = a*term;
			for (i = 0; i < 2; i++)
				pn[i+4] = b*pn[i+2] - an*pn[i];
			if (pn[5] != 0) {
				rn = pn[4]/pn[5];
				dif = fabs(gin-rn);
				if ((dif <= ACCURACY_TOL) && (dif <= ACCURACY_TOL*rn))
					return 1 - factor*gin;
				gin = rn;
			}
			for (i = 0; i < 4; i++)
				pn[i]= pn[i+2];
			if (fabs(pn[4]) >= OVERFLOW_VAL) {
				for (i = 0; i < 4; i++)
					pn[i] /= OVERFLOW_VAL;
			}
		}
	}
	else {
		/* series expansion */
		gin = 1;
		term = 1;
		rn = p;
		while (term > ACCURACY_TOL) {
			rn++;
			term *= x/rn;
			gin += term;
		}
		gin *= factor/p;
	}
	return (gin);
}
long factorial(int n) {
	long f, i;
	assert(n <= 10);
	for (i = 2, f = 1; i <= (long)n; i++)
		f*= i;
	return f;
}
double LnGamma (double x)
{
/* returns ln(gamma(x)) for x>0, accurate to 10 decimal places.
   Stirling's formula is used for the central polynomial part of the procedure.
   Pike MC & Hill ID (1966) Algorithm 291: Logarithm of the gamma function.
   Communications of the Association for Computing Machinery, 9:684
*/
	double f = 0, fneg = 0, z, tmp1, tmp2;
	int nx =(int)x-1;
	if((double)nx == x && nx > 0 && nx < 10)
		return log((double)factorial(nx));
	assert(x>0);
	if (x < 7) {
		f = 1;
		z = x-1;
		while (++z < 7)
			f *= z;
		x = z;
		f = -log(f);
	}
	z = 1/(x*x);
	tmp1 = fneg+ f + (x-0.5)*log(x) - x + .918938533204673;
	tmp2 = (((-.000595238095238*z+.000793650793651)*z-.002777777777778)*z +.083333333333333);
	return tmp1 + (tmp2/x);
}
/* functions concerning the CDF and percentage points of the gamma and
   Chi2 distribution
*/
double PointNormal (double prob) {
/* returns z so that Prob{x < z}= prob where x ~ N(0,1) and (1e-12)< prob < 1-(1e-12)
   returns (-9999) if in error
   Odeh RE & Evans JO (1974) The percentage points of the normal distribution.
   Applied Statistics 22: 96-97 (AS70)
   Newer methods:
	 Wichura MJ (1988) Algorithm AS 241: the percentage points of the
	   normal distribution.	 37: 477-484.
	 Beasley JD & Springer SG  (1977).	Algorithm AS 111: the percentage
	   points of the normal distribution.  26: 118-121.
*/
	double a0 =-.322232431088, a1 =-1, a2 =-.342242088547, a3 =-.0204231210245;
	double a4 =-.453642210148e-4, b0 =.0993484626060, b1 =.588581570495;
	double b2 =.531103462366, b3 =.103537752850, b4 =.0038560700634;
	double y, z = 0, p = prob, p1;
	p1 = (p < 0.5 ? p : 1-p);
	if (p1 < 1e-20)
		z = 999;
	else {
		y = sqrt (log(1/(p1*p1)));
		z = y + ((((y*a4+a3)*y+a2)*y+a1)*y+a0) / ((((y*b4+b3)*y+b2)*y+b1)*y+b0);
	}
	return (p < 0.5 ? -z : z);
}
double PointChi2 (double prob, double v) {
	/* returns z so that Prob{x < z}= prob where x is Chi2 distributed with df = v
	returns -1 if in error.	  0.000002 < prob < 0.999998
	RATNEST FORTRAN by
	Best DJ & Roberts DE (1975) The percentage points of the
	Chi2 distribution.	Applied Statistics 24: 385-388.	 (AS91)
	Converted into C by Ziheng Yang, Oct. 1993.
	*/
	double e =.5e-6, aa =.6931471805, p = prob, g, TINY_PROB = 1e-6;
	double xx, c, ch, a = 0,q = 0,p1 = 0,p2 = 0,t = 0,x = 0,b = 0,s1,s2,s3,s4,s5,s6;
	int doL3;
	assert(v > 0);
	if (p < TINY_PROB)
		return 0.0;
	if (p > 1.0 - TINY_PROB)
		return 9999.0;
	xx = v/2;
	g = LnGamma(xx);
	c = xx-1;
	if (v < -1.24*log(p)){
		ch = pow((p*xx*exp(g+xx*aa)), 1/xx);
		if (ch-e < 0)
			return (ch);
	}
	else {
		doL3 = 1;
		if (v <= .32) {
			ch = 0.4;
			a = log(1-p);
			for (;;) {
				q = ch;
				p1 = 1+ch*(4.67+ch);
				p2 = ch*(6.73+ch*(6.66+ch));
				t = -0.5+(4.67+2*ch)/p1 - (6.73+ch*(13.32+3*ch))/p2;
				ch -= (1-exp(a+g+.5*ch+c*aa)*p2/p1)/t;
				if (fabs(q/ch-1)-.01 <= 0) {
					doL3 = 0;
					break;
				}
			}
		}
		if (doL3) {
			x = PointNormal (p);
			p1 = 0.222222/v;
			ch = v*pow((x*sqrt(p1)+1-p1), 3.0);
			if (ch>2.2*v+6)
				ch = -2*(log(1-p)-c*log(.5*ch)+g);
		}
	}
	do {
		q = ch;
		p1 = .5*ch;
		t = IncompleteGamma (p1, xx, g);
		assert(t >= 0);
		p2 = p-t;
		t = p2*exp(xx*aa+g+p1-c*log(ch));
		b = t/ch;
		a = 0.5*t-b*c;
		s1 = (210+a*(140+a*(105+a*(84+a*(70+60*a))))) / 420;
		s2 = (420+a*(735+a*(966+a*(1141+1278*a))))/2520;
		s3 = (210+a*(462+a*(707+932*a)))/2520;
		s4 = (252+a*(672+1182*a)+c*(294+a*(889+1740*a)))/5040;
		s5 = (84+264*a+c*(175+606*a))/2520;
		s6 = (120+c*(346+127*c))/5040;
		ch += t*(1+0.5*t*s1-b*c*(s1-b*(s2-b*(s3-b*(s4-b*(s5-b*s6))))));
	} while(fabs(q/ch-1) > e);
	return ch;
}
/*
	Fills the first `n_cat` elements of `rates` with the mean rate of the
	corresponding quantile of the gamma distrib with rates alpha and gamma.
	The first element corresponds to the quantile demarcated by inverse CDF of
		0 up to 1/n_cat.
	The second goes from 1/n_cat up to 2/n_cat
	etc.
   Ziheng's comment on	DiscreteMean:
   discretization of gamma distribution with equal proportions in each
   category.
   MTH comment:
   I reworked this code from Ziheng's DiscreteMean to remove the need for a freqK
   array (which always ended up as n elements with values 1/n, but was used
   as scratch space internally).
*/
void DiscreteGammaMean (double rates[], const double alpha, const double beta, const unsigned n_cat)
{
	double lnga1;
	unsigned last_cat, i;
	double current;
	const double inv_two_beta = .5/beta;
	const double d_n_cat = (double)n_cat;
	const double cat_freq = 1.0/d_n_cat;
	double prop = cat_freq;
	double rate_sum = 0.0;
	double prev = 0.0;
	const double factor = (alpha*n_cat)/beta;
	assert(n_cat > 0);
	assert(rates != 0L);
	assert(alpha > 0.0);
	assert(beta > 0.0);
	if (n_cat < 2) {
		rates[0] = 1.0;
		return;
	}
	last_cat = n_cat - 1;
	lnga1 = LnGamma(alpha + 1.0);
	for (i = 0; i < last_cat; i++) {
		/*cutting points, Eq. 9 -- Ziheng's comment */
		current = inv_two_beta*PointChi2(prop, 2.0*alpha);
		current = IncompleteGamma(current*beta, alpha + 1, lnga1);
		rates[i] = (current - prev)*factor;
		rate_sum += rates[i];
		prev = current;
		prop +=	 cat_freq;
	}
	rates[last_cat] = (1.0 - prev)*factor;
	rate_sum += rates[last_cat];
	rate_sum /= d_n_cat;
	for (i = 0; i < n_cat; ++i)
		rates[i] /= rate_sum;
 }
/*
	Fills the first `n_cat` elements of `rates` with the mean rate of the
	corresponding quantile of the gamma distrib with rates alpha and gamma.
	The first element corresponds to the quantile demarcated by inverse CDF of
		0 up to 1/n_cat.
	The second goes from 1/n_cat up to 2/n_cat
	etc.
	Ziheng's comment on	DiscreteMean:
	discretization of gamma distribution with equal proportions in each
	category.
	MTH comment:
	I reworked this code from Ziheng's DiscreteMean to remove the need for a freqK
	array (which always ended up as n elements with values 1/n, but was used
	as scratch space internally).
*/
void DiscreteGammaMedian (double rates[], const double alpha, const double beta, unsigned n_cat)
{
	unsigned i;
	double total = 0.0;
	const double cat_freq = 1.0/((double)n_cat);
	double prop = cat_freq*0.5;
	double lnga1;
	assert(n_cat > 0);
	assert(rates != 0L);
	assert(alpha > 0.0);
	assert(beta > 0.0);
	if (n_cat < 2) {
		rates[0] = 1.0;
		return;
	}
	lnga1 = LnGamma(alpha + 1.0);
	for(i = 0; i < n_cat; ++i) {
		rates[i] = PointChi2(prop, 2.0*alpha);
		total += rates[i];
		prop +=	 cat_freq;
	}
	total /= (double)n_cat;
	for(i = 0; i < n_cat; ++i)
		rates[i] /= total;
}

void DiscreteGamma(double *f,double *r, double alpha, double beta, int ncat, int useMean) {
	const double freq = 1.0/((double)ncat);
	int i;
	for (i = 0; i < ncat; ++i)
		f[i] = freq;
	if (useMean)
		DiscreteGammaMean(r, alpha, beta, ncat);
	else
		DiscreteGammaMedian(r, alpha, beta, ncat);
}
/* end from Ziheng Yang's PAML */
/* end phylogenetic and numerical code */


/*---------------------------------------------------------------------------------
|	get_eigens
|
|	\Returns 0 on error, *is_complex will be 1 if the cause of the failure is
|		complex eigenvalues.
|	Taken from MrBayes GetEigens.
|	MrBayes seems to flag any complex eigensystems with a return code and abort
|		calculations.  Thus, I removed the complex ** arguments and code in the
|		isComplex branch.
|   Code taken from MrBayes 3.2 (see top of file)
---------------------------------------------------------------------------------*/
int get_eigens(
  unsigned dim,
  double **q,
  double *eigenValues,
  double *imEigenValues,
  double **eigenVectors,
  double **invEigenVectors,
  double **workMat, /*scratch  dim by dim matrix*/
  double *dWork, /*scratch of length dim*/
  int *iWork, /*scratch of length dim*/
  unsigned *is_complex) {
	int rc;
	*is_complex = 0;
	memset(dWork, 0, (size_t) (dim * sizeof(double)));
	memset(iWork, 0, (size_t) (dim * sizeof(int)));
	/* calculate eigenvalues and eigenvectors */
	rc = compute_eigen_system (dim, q, eigenValues, imEigenValues, eigenVectors, iWork, dWork, is_complex);
	if (rc == 0) {
		if (*is_complex == 1){
			PyErr_SetString(PyExc_ValueError,"Complex Eigenvalues found");
		}
		return 0;
	}
	CopyDoubleMatrices (dim, (const double **) eigenVectors, workMat);
	if (InvertMatrix (dim, workMat, dWork, iWork, invEigenVectors) == 0L)
		return 0;
	return 1;
}



/*---------------------------------------------------------------------------------
|   D_sign
|   This function is called from "Hqr2".
|   Code taken from MrBayes 3.2 (see top of file)
---------------------------------------------------------------------------------*/
double D_sign (double a, double b) {
	double		x;
	x = (a >= 0 ? a : -a);
	return (b >= 0 ? x : -x);
}
/*---------------------------------------------------------------------------------
|
|   ComplexDivision2
|
|   Returns the complex quotient of two complex numbers. It does not require that
|   the numbers be in a complex structure.
|
|   Code taken from MrBayes 3.2 (see top of file)
---------------------------------------------------------------------------------*/
void ComplexDivision2 (double ar, double ai, double br, double bi, double *cr, double *ci)
{
	double		s, ais, bis, ars, brs;
	s = fabs(br) + fabs(bi);
	ars = ar / s;
	ais = ai / s;
	brs = br / s;
	bis = bi / s;
	s = brs*brs + bis*bis;
	*cr = (ars*brs + ais*bis) / s;
	*ci = (ais*brs - ars*bis) / s;
}
/*---------------------------------------------------------------------------------
|
|   Exchange
|
|   Code taken from MrBayes 3.2 (see top of file)
---------------------------------------------------------------------------------*/
void Exchange (int j, int k, int l, int m, int n, double **a, double *scale)
{
	int			i;
	double		f;
	scale[m] = (double)j;
	if (j != m) {
		for (i = 0; i <= l; i++) {
			f = a[i][j];
			a[i][j] = a[i][m];
			a[i][m] = f;
		}
		for (i = k; i < n; i++) {
			f = a[j][i];
			a[j][i] = a[m][i];
			a[m][i] = f;
		}
	}
}
/*---------------------------------------------------------------------------------
|
|   ElmHes
|
|   Given a real general matrix, this subroutine
|   reduces a submatrix situated in rows and columns
|   low through high to upper Hessenberg form by
|   stabilized elementary similarity transformations.
|
|   On input:
|
|    * dim is the order of the matrix
|
|    * low and high are integers determined by the balancing
|      subroutine  balanc.  if  balanc  has not been used,
|      set low = 1, high = dim.
|
|    * a contains the input matrix.
|
|   On output:
|
|    * a contains the hessenberg matrix.  The multipliers
|      which were used in the reduction are stored in the
|      remaining triangle under the hessenberg matrix.
|
|    * interchanged contains information on the rows and columns
|      interchanged in the reduction.
|
|   Only elements low through high are used.
|
|   Code taken from MrBayes 3.2 (see top of file)
---------------------------------------------------------------------------------*/
void ElmHes (int dim, int low, int high, double **a, int *interchanged) {
	int			i, j, m, la, mm1, kp1, mp1;
	double		x, y;
	la = high - 1;
	kp1 = low + 1;
	if (la < kp1)
		return;
	for (m = kp1; m <= la; m++) {
		mm1 = m - 1;
		x = 0.0;
		i = m;
		for (j = m; j <= high; j++) {
			if (fabs(a[j][mm1]) > fabs(x)) {
				x = a[j][mm1];
				i = j;
			}
		}
		interchanged[m] = i;
		if (i != m)  {
			/* interchange rows and columns of a */
			for (j = mm1; j < dim; j++) {
				y = a[i][j];
				a[i][j] = a[m][j];
				a[m][j] = y;
			}
			for (j = 0; j <= high; j++) {
				y = a[j][i];
				a[j][i] = a[j][m];
				a[j][m] = y;
			}
		}
		if (fabs(x) > MB_ETA_TOL) {
			mp1 = m + 1;
			for (i = mp1; i <= high; i++) {
				y = a[i][mm1];
				if (fabs(y) > MB_ETA_TOL) {
					y /= x;
					a[i][mm1] = y;
					for (j = m; j < dim; j++)
						a[i][j] -= y * a[m][j];
					for (j = 0; j <= high; j++)
						a[j][m] += y * a[j][i];
				}
			}
		}
	}
}
/*---------------------------------------------------------------------------------
|
|   Balanc
|
|   This subroutine balances a real matrix and isolates
|   eigenvalues whenever possible.
|
|   On input:
|
|    * dim is the order of the matrix
|
|    * a contains the input matrix to be balanced
|
|   On output:
|
|    * a contains the balanced matrix.
|
|    * low and high are two integers such that a(i,j)
|      is equal to zero if
|         (1) i is greater than j and
|         (2) j = 1,...,low-1 or i = igh+1,...,n.
|
|    * scale contains information determining the
|      permutations and scaling factors used.
|
|   Suppose that the principal submatrix in rows pLow through pHigh
|   has been balanced, that p(j) denotes the index interchanged
|   with j during the permutation step, and that the elements
|   of the diagonal matrix used are denoted by d(i,j). Then
|      scale(j) = p(j),    for j = 1,...,pLow-1
|               = d(j,j),      j = pLow,...,pHigh
|               = p(j)         j = pHigh+1,...,dim.
|   The order in which the interchanges are made is dim to pHigh+1,
|   then 1 to pLow-1.
|
|   Note that 1 is returned for pHigh if pHigh is zero formally.
|
|   The algol procedure exc contained in balance appears in
|   balanc in line.  (Note that the algol roles of identifiers
|   k,l have been reversed.)
|
|   This routine is a translation of the Algol procedure from
|   Handbook for Automatic Computation, vol. II, Linear Algebra,
|   by Wilkinson and Reinsch, Springer-Verlag.
|
|   This function was converted from FORTRAN by D. L. Swofford.
|
|   Code taken from MrBayes 3.2 (see top of file)
---------------------------------------------------------------------------------*/
void Balanc (int dim, double **a, int *low, int *high, double *scale) {
	int			i, j, k, l, m, noconv;
	double		c, f, g, r, s, b2;
	b2 = FLT_RADIX * FLT_RADIX;
	k = 0;
	l = dim - 1;
	for (j = l; j>= 0; j--) {
		for (i = 0; i <= l; i++) {
			if ((i != j) && (fabs(a[j][i]) > MB_ETA_TOL)) {
				goto next_j1;
			}
		}
		/* bug that DLS caught */
		/*m = l;
		Exchange(j, k, l, m, dim, a, scale);
		if (l < 0)
			goto leave;
		else
			j = --l;*/
		m = l;
		Exchange(j, k, l, m, dim, a, scale);
		if (--l < 0)
			goto leave;
		next_j1:
			;
	}
	for (j = k; j <= l; j++) {
		for (i = k; i <= l; i++) {
			if ((i != j) && (fabs(a[i][j]) > MB_ETA_TOL)) {
				goto next_j;
			}
		}
		m = k;
		Exchange(j, k, l, m, dim, a, scale);
		k++;
		next_j:
			;
	}
	for (i = k; i <= l; i++)
		scale[i] = 1.0;
	noconv = 1;
	while (noconv) {
		noconv = 0;
		for (i = k; i <= l; i++) {
			c = 0.0;
			r = 0.0;
			for (j = k; j <= l; j++) {
				if (j != i) {
					c += fabs(a[j][i]);
					r += fabs(a[i][j]);
				}
			}
			if (fabs(c) > MB_ETA_TOL && fabs(r) > MB_ETA_TOL) {
				g = r / FLT_RADIX;
				f = 1.0;
				s = c + r;
				while (c < g) {
					f *= FLT_RADIX;
					c *= b2;
				}
				g = r * FLT_RADIX;
				while (c >= g) {
					f /= FLT_RADIX;
					c /= b2;
				}
				if ((c + r) / f < s * .95) {
					g = 1.0 / f;
					scale[i] *= f;
					noconv = 1;
					for (j = k; j < dim; j++)
						a[i][j] *= g;
					for (j = 0; j <= l; j++)
						a[j][i] *= f;
				}
			}
		}
	}
	leave:
		*low = k;
		*high = l;
}
/*---------------------------------------------------------------------------------
|
|   ElTran
|
|   This subroutine accumulates the stabilized elementary
|   similarity transformations used in the reduction of a
|   real general matrix to upper Hessenberg form by ElmHes.
|
|   On input:
|
|    * dim is the order of the matrix.
|
|    * low and high are integers determined by the balancing
|      subroutine  balanc. If Balanc has not been used,
|      set low = 0, high = dim-1.
|
|    * a contains the multipliers which were used in the
|      reduction by  ElmHes in its lower triangle
|      below the subdiagonal.
|
|    * interchanged contains information on the rows and columns
|      interchanged in the reduction by ElmHes.
|      only elements low through high are used.
|
|   On output:
|
|    * z contains the transformation matrix produced in the
|      reduction by ElmHes.
|
|   This routine is a translation of the Algol procedure from
|   Handbook for Automatic Computation, vol. II, Linear Algebra,
|   by Wilkinson and Reinsch, Springer-Verlag.
|
|   Code taken from MrBayes 3.2 (see top of file)
---------------------------------------------------------------------------------*/
void ElTran (int dim, int low, int high, double **a, int *interchanged, double **z)
{
	int			i, j, mp;
	/* initialize z to identity matrix */
	for (j = 0; j < dim; j++)
		{
		for (i = 0; i < dim; i++)
			z[i][j] = 0.0;
		z[j][j] = 1.0;
		}
	for (mp = high-1; mp>= low+1; mp--) /* there were a number of additional    */
		{                            /* variables (kl, la, m, mm, mp1) that  */
		for (i = mp+1; i <= high; i++)   /* have been eliminated here simply by  */
			z[i][mp] = a[i][mp-1];   /* initializing variables appropriately */
		i = interchanged[mp];        /* in the loops                         */
		if (i != mp) /* change "==" to "!=" to eliminate a goto statement */
			{
			for (j = mp; j <= high; j++)
				{
				z[mp][j] = z[i][j];
				z[i][j] = 0.0;
				}
			z[i][mp] = 1.0;
			}
		}
}
/*---------------------------------------------------------------------------------
|
|   Hqr2
|
|   This subroutine finds the eigenvalues and eigenvectors
|   of a real upper Hessenberg matrix by the QR method. The
|   eigenvectors of a real general matrix can also be found
|   if ElmHes  and ElTran or OrtHes and OrTran have
|   been used to reduce this general matrix to Hessenberg form
|   and to accumulate the similarity transformations.
|
|   On input:
|
|    * dim is the order of the matrix.
|
|    * low and high are integers determined by the balancing
|      subroutine  balanc. If  balanc has not been used,
|      set low = 0, high = dim-1.
|
|    * h contains the upper hessenberg matrix. Information about
|      the transformations used in the reduction to Hessenberg
|      form by  ElmHes  or OrtHes, if performed, is stored
|      in the remaining triangle under the Hessenberg matrix.
|
|   On output:
|
|    * h has been destroyed.
|
|    * wr and wi contain the real and imaginary parts,
|      respectively, of the eigenvalues. The eigenvalues
|      are unordered except that complex conjugate pairs
|      of values appear consecutively with the eigenvalue
|      having the positive imaginary part first. If an
|      error exit is made, the eigenvalues should be correct
|      for indices j,...,dim-1.
|
|    * z contains the transformation matrix produced by ElTran
|      after the reduction by ElmHes, or by OrTran after the
|      reduction by OrtHes, if performed. If the eigenvectors
|      of the Hessenberg matrix are desired, z must contain the
|      identity matrix.
|
|   Calls ComplexDivision2 for complex division.
|
|   This function returns:
|      zero       for normal return,
|      j          if the limit of 30*n iterations is exhausted
|                 while the j-th eigenvalue is being sought.
|
|   This subroutine is a translation of the ALGOL procedure HQR2,
|   Num. Math. 14, 219,231(1970) by Martin, Peters, and Wilkinson.
|   Handbook for Automatic Computation, vol. II - Linear Algebra,
|   pp. 357-391 (1971).
|
|   Code taken from MrBayes 3.2 (see top of file)
---------------------------------------------------------------------------------*/
int Hqr2 (int dim, int low, int high, double **h, double *wr, double *wi, double **z)
{
	int			i, j, k, l, m, na, en, notlas, mp2, itn, its, enm2, twoRoots;
	double		norm, p = 0.0, q = 0.0, r = 0.0, s = 0.0, t, w = 0.0, x, y = 0.0, ra, sa, vi, vr, zz = 0.0, tst1, tst2;
	norm = 0.0;
	k = 0;  /* used for array indexing. FORTRAN version: k = 1 */
	/* store roots isolated by balance, and compute matrix norm */
	for (i = 0; i < dim; i++)
		{
		for (j = k; j < dim; j++)
			norm += fabs(h[i][j]);
		k = i;
		if ((i < low) || (i > high))
			{
			wr[i] = h[i][i];
			wi[i] = 0.0;
			}
		}
	en = high;
	t = 0.0;
	itn = dim * 30;
	/* search for next eigenvalues */
	while (en >= low) /* changed from an "if(en < lo)" to eliminate a goto statement */
		{
		its = 0;
		na = en - 1;
		enm2 = na - 1;
		twoRoots = 0;
		for (;;)
			{
			for (l = en; l>low; l--) /* changed indexing, got rid of lo, ll */
				{
				s = fabs(h[l-1][l-1]) + fabs(h[l][l]);
				if (fabs(s) <= MB_ETA_TOL) /* == 0.0 */
					s = norm;
				tst1 = s;
				tst2 = tst1 + fabs(h[l][l-1]);
				if (fabs(tst2 - tst1) < MB_ETA_TOL) /* tst2 == tst1 */
					break; /* changed to break to remove a goto statement */
				}
			/* form shift */
			x = h[en][en];
			if (l == en) /* changed to break to remove a goto statement */
				break;
			y = h[na][na];
			w = h[en][na] * h[na][en];
			if (l == na)         /* used to return to other parts of the code */
				{
				twoRoots = 1;
				break;
				}
			if (itn == 0)
				return (en);
			/* form exceptional shift */
			if ((its == 10) || (its == 20)) /* changed to remove a goto statement */
				{
				t += x;
				for (i = low; i <= en; i++)
					h[i][i] -= x;
				s = fabs(h[en][na]) + fabs(h[na][enm2]);
				x = 0.75 * s;
				y = x;
				w = -0.4375 * s * s;
				}
			its++;
			itn--;
			/* look for two consecutive small sub-diagonal elements */
			for (m = enm2; m>= l; m--)
				{
				/* removed m = enm2 + l - mm and above loop to remove variables */
				zz = h[m][m];
				r = x - zz;
				s = y - zz;
				p = (r * s - w) / h[m+1][m] + h[m][m+1];
				q = h[m+1][m+1] - zz - r - s;
				r = h[m+2][m+1];
				s = fabs(p) + fabs(q) + fabs(r);
				p /= s;
				q /= s;
				r /= s;
				if (m == l)
					break; /* changed to break to remove a goto statement */
				tst1 = fabs(p) * (fabs(h[m-1][m-1]) + fabs(zz) + fabs(h[m+1][m+1]));
				tst2 = tst1 + fabs(h[m][m-1]) * (fabs(q) + fabs(r));
				if (fabs(tst2 - tst1) < MB_ETA_TOL) /* tst2 == tst1 */
					break; /* changed to break to remove a goto statement */
				}
			mp2 = m + 2;
			for (i = mp2; i <= en; i++)
				{
				h[i][i-2] = 0.0;
				if (i != mp2) /* changed "==" to "!=" to remove a goto statement */
					h[i][i-3] = 0.0;
				}
			/* double QR step involving rows l to en and columns m to en */
			for (k = m; k <= na; k++)
				{
				notlas = (k != na);
				if (k != m) /* changed "==" to "!=" to remove a goto statement */
					{
					p = h[k][k-1];
					q = h[k+1][k-1];
					r = 0.0;
					if (notlas)
						r = h[k+2][k-1];
					x = fabs(p) + fabs(q) + fabs(r);
					if (x < MB_ETA_TOL) /* == 0.0 */
						continue; /* changed to continue remove a goto statement */
					p /= x;
					q /= x;
					r /= x;
					}
		        /*
		        s = sqrt(p*p+q*q+r*r);
		        sgn = (p < 0)?-1:(p>0);
		        s = sgn*sqrt(p*p+q*q+r*r);
		        */
				s = D_sign(sqrt(p*p + q*q + r*r), p);
				if (k != m) /* changed "==" to "!=" to remove a goto statement */
					h[k][k-1] = -s * x;
				else if (l != m) /* else if gets rid of another goto statement */
					h[k][k-1] = -h[k][k-1];
				p += s;
				x = p / s;
				y = q / s;
				zz = r / s;
				q /= p;
				r /= p;
				if (!notlas) /* changed to !notlas to remove goto statement (see **) */
					{
					/* row modification */
					for (j = k; j < dim; j++)
						{
						p = h[k][j] + q * h[k+1][j];
						h[k][j] -= p * x;
						h[k+1][j] -= p * y;
						}
				    j = k + 3;
				    if (en < j)
				    	j = en;
					/* column modification */
					for (i = 0; i <= j; i++)
						{
						p = x * h[i][k] + y * h[i][k+1];
						h[i][k] -= p;
						h[i][k+1] -= p * q;
						}
					/* accumulate transformations */
					for (i = low; i <= high; i++)
						{
						p = x * z[i][k] + y * z[i][k+1];
						z[i][k] -= p;
						z[i][k+1] -= p * q;
						}
					}
				else /* (**) also put in else */
					{
					/* row modification */
					for (j = k; j < dim; j++)
						{
						p = h[k][j] + q * h[k+1][j] + r * h[k+2][j];
						h[k][j] -= p * x;
						h[k+1][j] -= p * y;
						h[k+2][j] -= p * zz;
						}
					j = k + 3;
					if (en < j)
						j = en;
					/* column modification */
					for (i = 0; i <= j; i++)
						{
						p = x * h[i][k] + y * h[i][k+1] + zz * h[i][k+2];
						h[i][k] -= p;
						h[i][k+1] -= p * q;
						h[i][k+2] -= p * r;
						}
					/* accumulate transformations */
					for (i = low; i <= high; i++)
						{
						p = x * z[i][k] + y * z[i][k+1] + zz * z[i][k+2];
						z[i][k] -= p;
						z[i][k+1] -= p * q;
						z[i][k+2] -= p * r;
						}
					}
				}
			}
		if (twoRoots)
			{
			/* two roots found */
			p = (y - x) / 2.0;
			q = p * p + w;
			zz = sqrt(fabs(q));
			h[en][en] = x + t;
			x = h[en][en];
			h[na][na] = y + t;
			if (q >= -1e-12) /* change "<" to ">=", and also change "0.0" to */
				{            /* a small number (Swofford's change)           */
				/* real pair */
				zz = p + D_sign(zz, p);
				wr[na] = x + zz;
				wr[en] = wr[na];
				if (fabs(zz) > MB_ETA_TOL) /* != 0.0 */
					wr[en] = x - w/zz;
				wi[na] = 0.0;
				wi[en] = 0.0;
				x = h[en][na];
				s = fabs(x) + fabs(zz);
				p = x / s;
				q = zz / s;
				r = sqrt(p*p + q*q);
				p /= r;
				q /= r;
				/* row modification */
				for (j = na; j < dim; j++)
					{
					zz = h[na][j];
					h[na][j] = q * zz + p * h[en][j];
					h[en][j] = q * h[en][j] - p * zz;
					}
				/* column modification */
				for (i = 0; i <= en; i++)
					{
					zz = h[i][na];
					h[i][na] = q * zz + p * h[i][en];
					h[i][en] = q * h[i][en] - p * zz;
					}
				/* accumulate transformations */
				for (i = low; i <= high; i++)
					{
					zz = z[i][na];
					z[i][na] = q * zz + p * z[i][en];
					z[i][en] = q * z[i][en] - p * zz;
					}
				}
			else
				{
				/* complex pair */
				wr[na] = x + p;
				wr[en] = x + p;
				wi[na] = zz;
				wi[en] = -zz;
				}
			en = enm2;
			}
		else
			{
			/* one root found */
			h[en][en] = x + t;
			wr[en] = h[en][en];
			wi[en] = 0.0;
			en = na;
			}
		}
	if (fabs(norm) < MB_ETA_TOL) /* == 0.0 */
		return (0); /* was a goto end of function */
	for (en = dim-1; en>= 0; en--)
		{
		/*en = n - nn - 1; and change for loop */
		p = wr[en];
		q = wi[en];
		na = en - 1;
		if (q < -1e-12)
			{
			/* last vector component chosen imaginary so that eigenvector
			   matrix is triangular */
			m = na;
			if (fabs(h[en][na]) > fabs(h[na][en]))
				{
				h[na][na] = q / h[en][na];
				h[na][en] = -(h[en][en] - p) / h[en][na];
				}
			else
				ComplexDivision2 (0.0, -h[na][en], h[na][na] - p, q, &h[na][na], &h[na][en]);
			h[en][na] = 0.0;
			h[en][en] = 1.0;
			enm2 = na - 1;
			if (enm2 >= 0) /* changed direction to remove goto statement */
				{
				for (i = enm2; i>= 0; i--)
					{
					w = h[i][i] - p;
					ra = 0.0;
					sa = 0.0;
					for (j = m; j <= en; j++)
						{
						ra += h[i][j] * h[j][na];
						sa += h[i][j] * h[j][en];
						}
					if (wi[i] < 0.0) /* changed direction to remove goto statement */
						{
						zz = w;
						r = ra;
						s = sa;
						}
					else
						{
						m = i;
						if (fabs(wi[i]) < MB_ETA_TOL) /* == 0.0 */ /* changed direction to remove goto statement */
							ComplexDivision2 (-ra, -sa, w, q, &h[i][na], &h[i][en]);
						else
							{
							/* solve complex equations */
							x = h[i][i+1];
							y = h[i+1][i];
							vr = (wr[i] - p) * (wr[i] - p) + wi[i] * wi[i] - q * q;
							vi = (wr[i] - p) * 2.0 * q;
							if ((fabs(vr) < MB_ETA_TOL) && (fabs(vi) < MB_ETA_TOL))
								{
								tst1 = norm * (fabs(w) + fabs(q) + fabs(x) + fabs(y) + fabs(zz));
								vr = tst1;
								do	{
									vr *= .01;
									tst2 = tst1 + vr;
									}
									while (tst2 > tst1); /* made into a do/while loop */
								}
							ComplexDivision2 (x * r - zz * ra + q * sa, x * s - zz * sa - q * ra, vr, vi, &h[i][na], &h[i][en]);
							if (fabs(x) > fabs(zz) + fabs(q)) /* changed direction to remove goto statement */
								{
								h[i+1][na] = (-ra - w * h[i][na] + q * h[i][en]) / x;
								h[i+1][en] = (-sa - w * h[i][en] - q * h[i][na]) / x;
								}
							else
								ComplexDivision2 (-r - y * h[i][na], -s - y * h[i][en], zz, q, &h[i+1][na], &h[i+1][en]);
							}
						/* overflow control */
						tst1 = fabs(h[i][na]);
						tst2 = fabs(h[i][en]);
						t = (tst1 > tst2 ? tst1 : tst2);
						if (t > MB_ETA_TOL) /* t != 0.0 */
							{
							tst1 = t;
							tst2 = tst1 + 1.0 / tst1;
							if (tst2 <= tst1)
								{
								for (j = i; j <= en; j++)
									{
									h[j][na] /= t;
									h[j][en] /= t;
									}
								}
							}
						}
					}
				}
			}
		else if (fabs(q)< MB_ETA_TOL)
			{
			/* real vector */
			m = en;
			h[en][en] = 1.0;
			if (na >= 0)
				{
				for (i = na; i>= 0; i--)
					{
					w = h[i][i] - p;
					r = 0.0;
					for (j = m; j <= en; j++)
						r += h[i][j] * h[j][en];
					if (wi[i] < 0.0) /* changed direction to remove goto statement */
						{
						zz = w;
						s = r;
						continue;  /* changed to continue to remove goto statement */
						}
					else
						{
						m = i;
						if (fabs(wi[i])< MB_ETA_TOL) /* changed to remove goto statement */
							{
							t = w;
							if (fabs(t)< MB_ETA_TOL)  /* changed to remove goto statement */
								{
								tst1 = norm;
								t = tst1;
								do	{
									t *= .01;
									tst2 = norm + t;
									}
									while (tst2 > tst1);
								}
							h[i][en] = -r / t;
							}
						else
							{
							/* solve real equations */
							x = h[i][i+1];
							y = h[i+1][i];
							q = (wr[i] - p) * (wr[i] - p) + wi[i] * wi[i];
							t = (x * s - zz * r) / q;
							h[i][en] = t;
							if (fabs(x) > fabs(zz))  /* changed direction to remove goto statement */
								h[i+1][en] = (-r - w * t) / x;
							else
								h[i+1][en] = (-s - y * t) / zz;
							}
						/* overflow control */
						t = fabs(h[i][en]);
						if (t > MB_ETA_TOL)
							{
							tst1 = t;
							tst2 = tst1 + 1. / tst1;
							if (tst2 <= tst1)
								{
								for (j = i; j <= en; j++)
									h[j][en] /= t;
								}
							}
						}
					}
				}
			}
		}
	for (i = 0; i < dim; i++)
		{
		if ((i < low) || (i > high)) /* changed to rid goto statement */
			{
			for (j = i; j < dim; j++)
				z[i][j] = h[i][j];
			}
		}
	/* multiply by transformation matrix to give vectors of original
	   full matrix */
	for (j = dim-1; j>= low; j--)
		{
		m = (j < high ? j : high);
		for (i = low; i <= high; i++)
			{
			zz = 0.0;
			for (k = low; k <= m; k++)
				zz += z[i][k] * h[k][j];
			z[i][j] = zz;
			}
		}
	return (0);
}
/*---------------------------------------------------------------------------------
|
|   BalBak
|
|   This subroutine forms the eigenvectors of a real general
|   matrix by back transforming those of the corresponding
|   balanced matrix determined by  balance.
|
|   On input:
|
|    * dim is the order of the matrix
|
|    * low and high are integers determined by  balance
|
|    * scale contains information determining the permutations
|      and scaling factors used by balance
|
|    * m is the number of columns of z to be back transformed
|
|    * z contains the real and imaginary parts of the eigen-
|      vectors to be back transformed in its first m columns
|
|   On output:
|
|    * z contains the real and imaginary parts of the
|      transformed eigenvectors in its first m columns
|
|   This routine is a translation of the Algol procedure from
|   Handbook for Automatic Computation, vol. II, Linear Algebra,
|   by Wilkinson and Reinsch, Springer-Verlag.
|
|   Code taken from MrBayes 3.2 (see top of file)
---------------------------------------------------------------------------------*/
void BalBak (int dim, int low, int high, double *scale, int m, double **z)
{
	int			i, j, k, ii;
	double		s;
	if (m != 0) /* change "==" to "!=" to eliminate a goto statement */
		{
		if (high != low) /* change "==" to "!=" to eliminate a goto statement */
			{
			for (i = low; i <= high; i++)
				{
				s = scale[i];
				for (j = 0; j < m; j++)
					z[i][j] *= s;
				}
			}
		for (ii = 0; ii < dim; ii++)
			{
			i = ii;
			if ((i < low) || (i > high)) /* was (i >= lo) && (i <= hi) but this */
				{                        /* eliminates a goto statement        */
				if (i < low)
					i = low - ii;
				k = (int)scale[i];
				if (k != i) /* change "==" to "!=" to eliminate a goto statement */
					{
					for (j = 0; j < m; j++)
						{
						s = z[i][j];
						z[i][j] = z[k][j];
						z[k][j] = s;
						}
					}
				}
			}
		}
}
/*---------------------------------------------------------------------------------
|
|   eigens_for_real_matrix
|
|   The matrix of interest is a. The ouptut is the real and imaginary parts of the
|   eigenvalues (wr and wi). z contains the real and imaginary parts of the
|   eigenvectors. iv2 and fv1 are working vectors.
|
|	Returns 0 on failure, 1 for success.
|
|   Code taken from MrBayes 3.2 (see top of file)
---------------------------------------------------------------------------------*/
int eigens_for_real_matrix (int dim, double **a, double *wr, double *wi, double **z, int *iv1, double *fv1)
{
	static int	is1, is2;
	int			ierr;
	Balanc (dim, a, &is1, &is2, fv1);
	ElmHes (dim, is1, is2, a, iv1);
	ElTran (dim, is1, is2, a, iv1, z);
	ierr = Hqr2 (dim, is1, is2, a, wr, wi, z);
	if (ierr == 0)
		BalBak (dim, is1, is2, fv1, dim, z);
	else {
		PyErr_SetString(PyExc_RuntimeError,"Error in Hqr2: max iterations exceeded.");
		return 0;
	}
	return 1;
}

/*---------------------------------------------------------------------------------
|
|   compute_eigen_system
|
|   Calculates the eigenvalues, eigenvectors, and the inverse of the eigenvectors
|   for a matrix of real numbers.
|
|	Returns 0 for an error or complex eigen values.  If the error is ONLY the
| 		presence of complex eigenvalues then, *is_complex = 1.
|
|   Code taken from MrBayes 3.2 (see top of file)
---------------------------------------------------------------------------------*/
static int compute_eigen_system (int dim, double **a, double *v, double *vi, double **u, int *iwork, double *dwork, unsigned *is_complex) {
	int i;
	*is_complex = 0;
	if (0 == eigens_for_real_matrix (dim, a, v, vi, u, iwork, dwork))
		return 0;
	for (i = 0; i < dim; i++) {
		if (fabs(vi[i]) > MB_ETA_TOL) { /* != 0.0 */
			*is_complex = 1;
			return 0;
		}
	}
	return 1;
}


/*---------------------------------------------------------------------------------
|   CopyDoubleMatrices
|   Copies the contents of one matrix of doubles to another matrix.
|   Code taken from MrBayes 3.2 (see top of file)
---------------------------------------------------------------------------------*/
static void CopyDoubleMatrices (int dim, const double **from, double **to) {

	int			i, j;
	for (i=0; i<dim; i++) {
		for (j=0; j<dim; j++) {
			to[i][j] = from[i][j];
		}
	}
}
/*---------------------------------------------------------------------------------
|
|   InvertMatrix
|
|   Calculates aInv = a^{-1} using LU-decomposition. The input matrix a is
|   destroyed in the process. The program returns an error if the matrix is
|   singular. dWork and iWork are work vectors.
|
|   Code taken from MrBayes 3.2 (see top of file)
---------------------------------------------------------------------------------*/
static int InvertMatrix (
	int dim,
 	double **a, /*dim by dim - valid on input, but overwritten */
 	double *dWork, /* scratch array, len dim */
 	int *iWork, /* scratch array, len dim */
 	double **aInv/*dim by dim - valid on completion unless 0 is returned */
 	) {
	unsigned	i, j;
	if (LUDecompose (dim, a, dWork, iWork) == 0L)
		return 0;
	for (j = 0; j < dim; j++)
		{
		for (i = 0; i < dim; i++)
			dWork[i] = 0.0;
		dWork[j] = 1.0;
		LUBackSubstitution (dim, a, iWork, dWork);
		for (i = 0; i < dim; i++)
			aInv[i][j] = dWork[i];
		}
	return 1;
}


/*---------------------------------------------------------------------------------
|
|   LUBackSubstitution
|
|   Back substitute into an LU-decomposed matrix.
|
|   Code taken from MrBayes 3.2 (see top of file)
---------------------------------------------------------------------------------*/
static void LUBackSubstitution (int dim, double **a, int *iWork, double *b) {
	int			i, ip, j, ii = -1;
	double		sum;
	for (i = 0; i < dim; i++) {
		ip = iWork[i];
		sum = b[ip];
		b[ip] = b[i];
		if (ii >= 0) {
			for (j = ii; j <= i-1; j++)
				sum -= a[i][j] * b[j];
		}
		else if (fabs(sum) > MB_ETA_TOL)
			ii = i;
		b[i] = sum;
	}
	for (i = dim-1; i>= 0; i--) {
		sum = b[i];
		for (j = i+1; j < dim; j++)
			sum -= a[i][j] * b[j];
		b[i] = sum / a[i][i];
	}
}
/*---------------------------------------------------------------------------------
|
|   LUDecompose
|
|   Calculate the LU-decomposition of the matrix a. The matrix a is replaced.
|
|   Code taken from MrBayes 3.2 (see top of file)
---------------------------------------------------------------------------------*/
static int LUDecompose (int dim, double **a, double *dWork, int *iWork) {
	int			i, imax = 0, j, k;
	double		big, dum, sum, temp, d;
	d = 1.0;
	for (i = 0; i < dim; i++) {
		big = 0.0;
		for (j = 0; j < dim; j++) {
			temp = fabs(a[i][j]);
			if (temp > big)
				big = temp;
			}
		if (fabs(big) < MB_ETA_TOL) {
			PyErr_SetString(PyExc_RuntimeError,"Error in LUDecompose: largest abs value of the elements in a row is too small.");
			return 0;
		}
		dWork[i] = 1.0 / big;
	}
	for (j = 0; j < dim; j++) {
		for (i = 0; i < j; i++) {
			sum = a[i][j];
			for (k = 0; k < i; k++)
				sum -= a[i][k] * a[k][j];
			a[i][j] = sum;
		}
		big = 0.0;
		for (i = j; i < dim; i++) {
			sum = a[i][j];
			for (k = 0; k < j; k++)
				sum -= a[i][k] * a[k][j];
			a[i][j] = sum;
			dum = dWork[i] * fabs(sum);
			if (dum >= big) {
				big = dum;
				imax = i;
			}
		}
		if (j != imax) {
			for (k = 0; k < dim; k++) {
				dum = a[imax][k];
				a[imax][k] = a[j][k];
				a[j][k] = dum;
			}
			d = -d;
			dWork[imax] = dWork[j];
		}
		iWork[j] = imax;
		if (fabs(a[j][j]) < MB_ETA_TOL)
			a[j][j] = MB_TINY_TOL;
		if (j != dim - 1) {
			dum = 1.0 / a[j][j];
			for (i = j + 1; i < dim; i++)
				a[i][j] *= dum;
		}
	}
	return 1;
}

