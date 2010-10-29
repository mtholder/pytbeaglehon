#include "pytbeaglehon/ccore/asrv.h"
#include "pytbeaglehon/ccore/calc_instance.h"
#include "pytbeaglehon/ccore/discrete_state_model.h"
#include <stdlib.h>
#include <stdio.h>
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
    fprintf(stderr, "%d failures out of %d tests.\n", numErrors, numErrors + numPasses);
    return numErrors;
}
