#! /usr/bin/env python
# Copyright (c) 2005-7 by Mark T. Holder
# (see bottom of file)

"unit tests of among-site rate heterogeneity code"
import unittest
from pytbeaglehon.tests.util import assert_list_eq
# pylint: disable-msg=C0111,W0401,W0611,W0212
from pytbeaglehon.parameter import Parameter, FloatParameter, MutableFloatParameter, ProbabilityVectorParameter
from pytbeaglehon import get_logger

_LOG = get_logger(__name__)

class ParameterTest(unittest.TestCase):
    def test_simple(self):
        fp = FloatParameter(1.0)
        self.assertRaises(AttributeError, fp.set_value, 3)
        mfp = MutableFloatParameter(1.0)
        self.assertEqual(mfp.value, 1.0)
        mfp.value = 2.5
        self.assertEqual(mfp.value, 2.5)
        self.assertRaises(ValueError, mfp.set_value, "bogus")
        
# pylint: disable-msg=C0103
def getTestSuite():
    """Alias to the additional_tests().  This is unittest-style.
    additional_tests is used by setuptools.
    """
    return additional_tests()

if __name__ == "__main__":
    unittest.main()

##############################################################################
##  pytbeaglehon phylogenetic likelihood caluclations using beaglelib.
##
##  Copyright 2010 Mark T. Holder
##  All rights reserved.
##
##  See "LICENSE.txt" for terms and conditions of usage.
##
##############################################################################
