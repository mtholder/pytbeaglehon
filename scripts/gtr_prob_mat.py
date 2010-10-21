#!/usr/bin/env python
from pytbeaglehon.disc_state_cont_time_model import RevDiscStateContTimeModel


def demand_positive_float(arg, param_name):
    try:
        v = float(arg)
        if v <= 0.0:
            sys.exit("Expecting %s to be greater than 0, but found %s." % (param_name, arg))
        return v
    except:
        sys.exit("Expecting a number for %s, but found %s." % (param_name, arg))

help_message = """gtr_prob_mat prints out the transition probablity matrix for the GTR model.

Options:
-h                                This help message
-f <freq_A> <freq_C> <freq_G>     Base frequencies
-r <rAC> <rAG> <rAT> <rCG> <rCT>  The relative rate matrix
-e <edge length>                  The edge length
"""
import sys
from itertools import izip
freq = [.25]*4
r_mat = [[1.0, 1.0, 1.0], [1.0, 1.0], [1.0]]
edge_len = 0.01

reading_freq, reading_rmat, reading_edge_len = False, False, False
read_freq, read_rmat, read_edge_len = False, False, False
for arg in sys.argv[1:]:
    if reading_freq:
        curr_list.append(demand_positive_float(arg, 'base frequency'))
        if len(curr_list) > 2:
            s = sum(curr_list)
            if s > 1.0:
                sys.exit("Expecting the frequency of A, C, and G to add up to less than 1")
            freq = curr_list + [1.0 - s]
            curr_list = []
            reading_freq = False
            read_freq = True
    elif reading_rmat:
        curr_list.append(demand_positive_float(arg, 'relative rate'))
        if len(curr_list) > 4:
            r_mat = [curr_list[:3], curr_list[3:], [1.0]]
            curr_list = []
            reading_rmat = False
            read_rmat = True
    elif reading_edge_len:
        edge_len = demand_positive_float(arg, 'edge length')
        read_edge_len = True
        reading_edge_len = False
    else:
        if arg == '-h':
            sys.stderr.write(help_message)
            sys.exit(0)
        if arg == '-f':
            if read_freq:
                sys.exit("Expected only one -f per invocation")
            curr_list = []
            reading_freq = True
        elif arg == '-r':
            if read_rmat:
                sys.exit("Expected only one -r per invocation")
            curr_list = []
            reading_rmat = True
        elif arg == '-e':
            if read_edge_len:
                sys.exit("Expected only one -e per invocation")
            curr_list = []
            reading_edge_len = True

if reading_freq:
    sys.exit("Expected 3 base frequencies")
elif reading_rmat:
    sys.exit("Expected 5 relative rate parameters")
elif reading_edge_len:
    sys.exit("Expected an edge length parameter")

# Arguments have been parsed. Do the calculation...

model = RevDiscStateContTimeModel(state_freq=freq, r_upper=r_mat)
mat = model.calc_prob_matrices(edge_len)[0] # returns a list of matrices, but we are not 
                                            #   using rate het, so we just want
                                            #   the first.


# format the output

max_len = 0
mat_str = [list(' ACGT')]
for from_base, mat_row in izip('ACGT', mat):
    mat_str_row = [from_base]
    for cell in mat_row:
        s = str(cell)
        if len(s) > max_len:
            max_len = len(s)
        mat_str_row.append(s)
    mat_str.append(mat_str_row)
format_str = "%%-%ds" % max_len
for row in mat_str:
    print "   ".join([format_str % i for i in row])
    
##############################################################################
##  pytbeaglehon phylogenetic likelihood caluclations using beaglelib.
##
##  Copyright 2010 Mark T. Holder
##  All rights reserved.
##
##  See "LICENSE.txt" for terms and conditions of usage.
##
##############################################################################
