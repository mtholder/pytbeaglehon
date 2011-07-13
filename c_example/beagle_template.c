#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <libhmsbeagle/beagle.h>
#include "pytbeaglehon/ccore/phylo_util.h"
#include "pytbeaglehon/ccore/asrv.h"
#include "pytbeaglehon/ccore/calc_instance.h"
#include "pytbeaglehon/ccore/discrete_state_model.h"
#include "pytbeaglehon/ccore/internal_like_calc_env.h"

void writeProbMat(int instance, int pmatind, double ** pmatBuffer, double edgeLen) {
    int j, k;
    beagleGetTransitionMatrix(instance, pmatind, pmatBuffer[0]);
    fprintf(stderr, "ProbMat #%d.  Edgelen = %lf", pmatind, edgeLen);
    for (k = 0; k < 4; ++k) {
        fprintf(stderr, "\n  ");
        for (j = 0; j < 4; ++j) {
            fprintf(stderr, "%lf, ", pmatBuffer[k][j]);
        }
    }
    fprintf(stderr, "\n");
}


void writePartial(int instance, int partind, double *partialBuffer, int nstates, int nchars) {
    int j, k;
    beagleGetPartials(instance, partind, BEAGLE_OP_NONE, partialBuffer);
    fprintf(stderr, "Partial #%d.", partind);
    for (k = 0; k < nchars; ++k) {
        fprintf(stderr, "\n  char %d => ", k);
        for (j = 0; j < nstates ; ++j) 
            fprintf(stderr, "%lf, ", *partialBuffer++);
    }
    fprintf(stderr, "\n");
}



int main(int argc, char * argv[]) {
    BeagleInstanceDetails beagleInstanceDetails ; 
    long resourcePref; 
    long resourceReq;
    int rc;
    int resourceIndex;
    fprintf(stderr, "Running BEAGLE_API test\n");

/*  Code for printing out a few probmats (hard codede in inds and numC)
    int i;
    int inds[] = {12, 0, 9, 11, 4, 2};
    int numC = 6;
    int pmatind;
    for (i = 0; i < numC; ++i) {
        pmatind = inds[i];
        writeProbMat(0, pmatind, probMatTemp, inst0edgeLenArray[i]);
    }
*/
