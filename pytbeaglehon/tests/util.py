#! /usr/bin/env python
'''Utility functions to make it easier to write tests'''
# Copyright (c) 2010 by Mark T. Holder
# (see bottom of file)
from itertools import izip
import unittest
from pytbeaglehon import get_logger
_LOG = get_logger(__name__)

def assert_list_of_mat_eq(self, returned, expected):
    for ret_mat, exp_mat in izip(returned, expected):
        assert_mat_eq(self, ret_mat, exp_mat)
    
def assert_mat_eq(self, returned, expected):
    "calls self.assertAlmostEqual for all elements of two matrices."
    _LOG.debug("returned = %s\nexpected = %s" %(str(returned), str(expected)))
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
