#include "pytbeaglehon/ccore/asrv.h"
#include "pytbeaglehon/ccore/calc_instance.h"
#include "pytbeaglehon/ccore/discrete_state_model.h"
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


void testCalcInstances(unsigned * numPasses, unsigned * numErrors) {
    char n[81];
    char description[81];
    long optsFlags, reqFlags, prefFlags, instanceHandle1, instanceHandle2;
    int rc;
    unsigned i, numModels;
    const DSCTModelObj ** modArray = 0L;
    if (gVerbose)
        fprintf(stderr, "Testing CalcInstances\n");
    i = getNumComputationalResources();
    fprintf(stderr, "%d computational resource(s) available\n", i);
    if (i == 0) {
        *numErrors += 1;
        fprintf(stderr, "%d computational resources available\n", i);
    }
    else
        *numPasses += 1;
    rc = getComputationalResourceDetails(0, n, description, &optsFlags, &reqFlags);
    if (rc == 0) {
        *numPasses += 1;
        if (gVerbose)
            fprintf(stderr, "First resource details:\n  name = %s\n  desc = %s\n  optFlags  = %ld\n  reqFlags = %ld\n", n, description, optsFlags, reqFlags);
    }
    else {
        *numErrors += 1;
        fprintf(stderr, "getComputationalResourceDetails returned %d\n", rc);
    }
    rc = getComputationalResourceDetails(i, n, description, &optsFlags, &reqFlags);
    if (rc != 0) {
        *numPasses += 1;
    }
    else {
        *numErrors += 1;
        fprintf(stderr, "getComputationalResourceDetails did not return an error code when resource %d was used as an arg.\n", i);
    }


    instanceHandle1 = createLikelihoodCalcInstance(4, /* leaves*/
            10, /* patterns */
            0L, /* pattern weights*/
            4, /* states */
            4, /* stateCodeArrays */
            10,  /* partials */
            1, /* distinct models (q-matrices) to allocated */
            0L, /* ASRVObj array */
            10, /* transition probability matrices */ 
            2, /* num of eigen solutions stored */ 
            1, /* num rescalers */
            0, /* the index of the computational resource to use */
            prefFlags,
            reqFlags); /* the beagle flags (see above) required of the computational resource */
    if (instanceHandle1 >= 0) {
        *numPasses += 1;
    }
    else {
        *numErrors += 1;
        fprintf(stderr, "createLikelihoodCalcInstance returned error code %d.\n", i);
    }
    instanceHandle2 = createLikelihoodCalcInstance(4, /* leaves*/
            10, /* patterns */
            0L, /* pattern weights*/
            4, /* states */
            4, /* stateCodeArrays */
            10,  /* partials */
            1, /* distinct models (q-matrices) to allocated */
            0L, /* ASRVObj array */
            10, /* transition probability matrices */ 
            2, /* num of eigen solutions stored */ 
            1, /* num rescalers */
            -1, /* -1 should mean "use any resoure" */
            prefFlags,
            reqFlags); /* the beagle flags (see above) required of the computational resource */
    if (instanceHandle2 >= 0) {
        *numPasses += 1;
    }
    else {
        *numErrors += 1;
        fprintf(stderr, "createLikelihoodCalcInstance returned error code %d.\n", i);
    }
    if (instanceHandle2 != instanceHandle1) {
        *numPasses += 1;
    }
    else {
        *numErrors += 1;
        fprintf(stderr, "duplicate instance handles returned by createLikelihoodCalcInstance.\n");
    }

    modArray = getModelList(instanceHandle1, &numModels);
    if (modArray == 0L || numModels != 1) {
            *numErrors += 1;
            fprintf(stderr, "could not obtain model list.\n");
    }
    else {
        *numPasses += 1;
    }

    
    
    if (instanceHandle1 >= 0) {
        if (freeLikeCalculatorInstance(instanceHandle1) < 0) {
            *numErrors += 1;
            fprintf(stderr, "freeLikeCalculatorInstance failed.\n");
        }   
        else
            *numPasses += 1;
    }
    if (instanceHandle2 >= 0) {
        if (freeLikeCalculatorInstance(instanceHandle2) < 0) {
            *numErrors += 1;
            fprintf(stderr, "freeLikeCalculatorInstance failed.\n");
        }   
        else
            *numPasses += 1;
    }
    if (instanceHandle1 >= 0) {
        if (freeLikeCalculatorInstance(instanceHandle1) >= 0) {
            *numErrors += 1;
            fprintf(stderr, "double call to freeLikeCalculatorInstance did not fail.\n");
        }
        else
            *numPasses += 1;
    }
    
}

void testASRV(unsigned * numPasses, unsigned * numErrors) {
    if (gVerbose)
        fprintf(stderr, "Testing ASRV\n");
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
    testCalcInstances(&numPasses, &numErrors);
    fprintf(stderr, "%d failures out of %d tests.\n", numErrors, numErrors + numPasses);
    return numErrors;
}
