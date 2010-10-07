#! /usr/bin/env python
'''Utility functions to make it easier to write tests'''
# Copyright (c) 2010 by Mark T. Holder
# (see bottom of file)
from itertools import izip
import unittest
def assert_mat_eq(self, returned, expected):
    "calls self.assertAlmostEqual for all elements of two matrices."
    for ret_row, exp_row in izip(returned, expected):
        assert_list_eq(self, ret_row, exp_row)

def assert_list_eq(self, returned, expected):
    "calls self.assertAlmostEqual for all elements of two lists."
    for ret_val, exp_val in izip(returned, expected):
        self.assertAlmostEqual(ret_val, exp_val, places=5)

def additional_tests():
    "returns all tests in this file as suite"
    return unittest.TestSuite([])


##############################################################################
##  pytbeaglehon phylogenetic likelihood caluclations using beaglelib.
##
##  Copyright 2010 Mark T. Holder
##  All rights reserved.
##
##  See "LICENSE.txt" for terms and conditions of usage.
##
##############################################################################
