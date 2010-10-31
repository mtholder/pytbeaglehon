#include <stdlib.h>
#include <stdio.h>
#include <libhmsbeagle/beagle.h>

int main(int argc, char * argv[]) {
    BeagleInstanceDetails beagleInstanceDetails ; 
    long resourcePref; 
    long resourceReq;
    int rc;
    int resourceIndex;
    fprintf(stderr, "Running BEAGLE_API test\n");
