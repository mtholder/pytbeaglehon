#! /usr/bin/env python
# Copyright (c) 2005-7 by Mark T. Holder
# (see bottom of file)

"unit tests of among-site rate heterogeneity code"
import unittest
from pytbeaglehon.tests.util import assert_list_eq
# pylint: disable-msg=C0111,W0401,W0611,W0212
from pytbeaglehon import LikeCalcEnvironment as LCE
from pytbeaglehon import get_logger

_LOG = get_logger(__name__)

class LCETest(unittest.TestCase):
    def get_inst(self):
        lc_env = LCE(num_leaves=6,
            num_patterns=10,
            #pattern_weight_list
            num_states=4,
            num_state_code_arrays=5,
            num_partials=11,
            num_model_matrices=2,
            #asrv_list
            num_prob_matrices=3,
            num_eigen_storage_structs=7,
            num_rescalings_multipliers=8
            )
        return lc_env
    def test_resource_info(self):
        t = LCE.query_comp_resource_info(0)
        _LOG.debug("info tuple = %s" % str(t))
        self.assertEqual(len(t), 4)
        lc_env = self.get_inst()
        self.assertRaises(ValueError, lc_env.get_comp_resource_info)
        lc_env.resource_index = 0
        t = lc_env.comp_resource_info
        self.assertEqual(len(t), 4)
        _LOG.debug("info tuple = %s" % str(t))
        self.assertRaises(IndexError, LCE.query_comp_resource_info, -1)

    def test_lazy_init(self):
        lc_env = LCE()
        self.assertFalse(lc_env._incarnated)

    def test_init(self):
        lc_env = self.get_inst()
        lc_env._do_beagle_init() # trigger initialization
        self.assertTrue(lc_env._incarnated)

    def test_num_resources(self):
        self.assertTrue(LCE.get_num_comp_resources() > 0)
        
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
