#include "pytbeaglehon/ccore/asrv.h"
#include "pytbeaglehon/ccore/calc_instance.h"
#include "pytbeaglehon/ccore/discrete_state_model.h"
#include "pytbeaglehon/ccore/internal_like_calc_env.h"
#include <stdlib.h>
#include <stdio.h>

int eigenIndex, numToCalc;
DSCTModelObj * dsct_model_obj = 0L;
double lnL;
int main(int argc, char * argv[]) {
    int rc;
