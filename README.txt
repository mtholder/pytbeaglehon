PYTbeagleHON is a python wrapper around the beaglelib library for likelihood
    computations in the domain of phylogenetics.
    
PYTbeagleHON is written by Mark T. Holder (see License), but it also provides
    extensions to the beagle API. These extensions include code taken from PAML
    (by Ziheng Yang) for calculations associated with gamma-distributed 
    rate-heterogeneity, and code taken from MrBayes (by John Huelsenbeck, 
    Fredrik Ronquist, and Paul van der Mark) and from David L. Swofford for 
    calculating eigensolutions.
    


The configure and Makefiles are for building a C library that wraps and augments
    beaglelib.
    

To build for python see the setup.py file (and its notes). Below is an example
    of a build script that demonstrates some of the environmental settings used
    by the build.





################################################################################
# Example full_rebuild.sh script (for mac build)
################################################################################
set -x
rm -rf ./build/lib.macosx-10.4-x86_64-2.6
rm -rf ./build/temp.macosx-10.4-x86_64-2.6
find pytbeaglehon -name "*.pyc" -exec rm {} \; 
rm ./pytbeaglehon/ccore/disc_state_cont_time_model.so



export BEAGLE_PREFIX=/home/joe.user/beagle-lib/buildcpu/installed
export NCL_PREFIX="/home/joe.user/ncl/branches/buildDynamic/installed"
export DYLD_LIBRARY_PATH="${DYLD_LIBRARY_PATH}:${BEAGLE_PREFIX}/lib:${NCL_PREFIX}/lib/ncl"
export PYTBEAGLEHON_LOGGING_LEVEL=debug
export CFLAGS="$CFLAGS -arch x86_64"
export PYTHONPATH=`pwd`


python setup.py build --use-ncl || exit
cp build/lib.macosx-10.4-x86_64-2.6/pytbeaglehon/ccore/disc_state_cont_time_model.so pytbeaglehon/ccore/




################################################################################
# Generating C code from python execution
################################################################################
Occassionally (typically for debugging purposes), it is useful to generate the 
C pytbeaglehon calls that the python code invokes (or the beagle calls). There is
some support for doing this.

If you use the invocation:
    $ python setup.py build --debug --trace-printing 
to build, then you can use the following shell script to create real_beagletrace.c
and real_pytbeaglehontrace.c in the c_example directory (at the time of writing
not all of the beagle calls are tracked, but the important pytbeaglehon calls are).

################################################################################
#!/bin/sh
set -x
if test -z $1
then
    echo "Expecting a path to an example python script" 
    exit 1
fi

cp c_example/beagle_template.c c_example/real_beagletrace.c
python "$1" 2>&1 | grep '/\* BEAGLE_API' >> c_example/real_beagletrace.c
echo 'return 0;}' >> c_example/real_beagletrace.c

cp c_example/template.c c_example/real_pytbeaglehontrace.c
python "$1" 2>&1 | grep '/\* c' >> c_example/real_pytbeaglehontrace.c
echo 'return 0;}' >> c_example/real_pytbeaglehontrace.c
################################################################################


