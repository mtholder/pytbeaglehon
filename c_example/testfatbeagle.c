#include "pytbeaglehon/ccore/asrv.h"
#include <stdlib.h>
#include <stdio.h>
const double TOL = 1e-5;

int gVerbose = 0;
int is_approx_equal(double x, double y) {
    double d = x - y;
    if (d < 0.0)
        return (d > -TOL ? 1 : 0);
    return (d < TOL ? 1 : 0);
}

int test_approx_equal(double x, double y) {
    if (is_approx_equal(x, y) == 1)
        return 1;
    fprintf(stderr, "%f != %f\n", x, y);
    return 0;
}

int array_test_eq(unsigned dim, const double * x , const double * y, unsigned * numPasses, unsigned * numErrors);


int array_test_eq(unsigned dim, const double * x , const double y [], unsigned * numPasses, unsigned * numErrors) {
    unsigned origNumErrors = *numErrors;
    unsigned i = 0;
    for (; i < dim; ++i) {
        if (test_approx_equal(x[i], y[i]) == 0)
            *numErrors += 1;
        else
            *numPasses += 1;
    }
    return (origNumErrors == *numErrors ? 1 : 0);
}

ASRVObj* asrv_obj_new(unsigned dim, int style, double param);
void asrv_obj_dtor(ASRVObj* asrh);
void internal_asrv_set_shape(ASRVObj *asrh, double val);


void testASRV(unsigned * numPasses, unsigned * numErrors) {
    if (gVerbose)
        fprintf(stderr, "Testing ASRV\n");
    unsigned i;
    unsigned ncat = 4;
    ASRVObj * asrv = asrv_obj_new(ncat, 1, 0.4);
    
    if (asrv == 0L) {
        *numErrors += 1;
            fprintf(stderr, "Could not allocate asrv object\n");
    }
    else {
        *numPasses += 1;
        const double freq[] = {0.25, 0.25, 0.25, 0.25};
        array_test_eq(4, asrv->freq, freq, numPasses, numErrors);
        const double rates[] = {0.01671441,  0.18175607, 0.73128067, 3.07024885};
        array_test_eq(4, asrv->rate, rates, numPasses, numErrors);

        internal_asrv_set_shape(asrv, 0.5);

        array_test_eq(4, asrv->freq, freq, numPasses, numErrors);
        const double rates2[] = {0.03338775, 0.25191592, 0.82026848, 2.89442785};
        array_test_eq(4, asrv->rate, rates2, numPasses, numErrors);

        asrv_obj_dtor(asrv);
    }
}
int main(int argc, char * argv[]) {
    unsigned numErrors = 0;
    unsigned numPasses = 0;
    unsigned i;
    for (i = 1; i < argc; ++i) {
        if (argv[i][0] == '-') {
            if (argv[i][1] == 'v')
                gVerbose = 1;
        }
    }
    testASRV(&numPasses, &numErrors);
    
    fprintf(stderr, "%d failures out of %d tests.\n", numErrors, numErrors + numPasses);
    return numErrors;
}
