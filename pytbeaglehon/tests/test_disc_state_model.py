#! /usr/bin/env python
# Copyright (c) 2005-7 by Mark T. Holder,  University of Kansas
# (see bottom of file)

"unit tests datatyp internals"
import unittest
from pytbeaglehon import get_logger
_LOG = get_logger(__name__)

from pytbeaglehon.tests.util import *
# pylint: disable-msg=C0111,W0401,W0611,W0212
from pytbeaglehon.parameter import MutableFloatParameter
from pytbeaglehon.disc_state_cont_time_model import JukesCantorModel, DNAType, RevDiscStateContTimeModel, Kimura2ParameterModel, HKY85Model

from pytbeaglehon.disc_state_cont_time_model import _r_upper_to_r_mat, _r_mat_to_r_upper

class ModelTest(unittest.TestCase):
    def test_gtr(self):
        kappa = 6.122449
        m = RevDiscStateContTimeModel(state_freq=(0.3, 0.25, 0.2, 0.25),
                                      r_upper=[[1.0, kappa,   1.0 ],
                                                    [  1.0, kappa ],
                                                           [  1.0 ]]
                                     )
        assert_list_of_mat_eq(self, m.calc_prob_matrices(0.01), [[[0.991443333333, 0.001247, 0.00606333333333, 0.001247, ],
                                            [0.0014964, 0.989928, 0.0009976, 0.007576],
                                            [0.009095, 0.001247, 0.988415, 0.001247],
                                            [0.0014964, 0.007576, 0.0009976, 0.989928]]])
    def test_modhky(self):
        m = HKY85Model(2.0, [.25, 0.25, 0.25, 0.25])
        m.kappa.is_mutable = True
        m.kappa = 6.122449
        m.state_freq.is_mutable = True
        m.state_freq = (0.3, 0.25, 0.2, 0.25)
        assert_list_of_mat_eq(self, m.calc_prob_matrices(0.01), [[[0.991443333333, 0.001247, 0.00606333333333, 0.001247, ],
                                            [0.0014964, 0.989928, 0.0009976, 0.007576],
                                            [0.009095, 0.001247, 0.988415, 0.001247],
                                            [0.0014964, 0.007576, 0.0009976, 0.989928]]])
    def test_change_kappa(self):
        mp = MutableFloatParameter(2.0)
        m = Kimura2ParameterModel(mp)
        nc, ti, tv = 0.951679099289, 0.0239356129609, 0.0121926438748
        assert_list_of_mat_eq(self, m.calc_prob_matrices(0.05), [[[nc, tv, ti, tv], [tv, nc, tv, ti], [ti, tv, nc, tv], [tv, ti, tv, nc]]])
        mp.value = 1.0
        nc, ti, tv = 0.951630238774, 0.0161232537421, 0.0161232537421
        assert_list_of_mat_eq(self, m.calc_prob_matrices(0.05), [[[nc, tv, ti, tv], [tv, nc, tv, ti], [ti, tv, nc, tv], [tv, ti, tv, nc]]])
        assert_list_of_mat_eq(self, m.calc_prob_matrices(0.05), [[[nc, tv, ti, tv], [tv, nc, tv, ti], [ti, tv, nc, tv], [tv, ti, tv, nc]]])
    def test_jc_probs(self):
        jc = JukesCantorModel()
        nc, c = 0.99006637135539677, 0.0033112095482010773
        assert_list_of_mat_eq(self, jc.calc_prob_matrices(0.01), [[[nc, c, c, c], [c, nc, c, c], [c, c, nc, c], [c, c, c, nc]]])
        nc, c = 1.0, 0.0
        assert_list_of_mat_eq(self, jc.calc_prob_matrices(0.0), [[[nc, c, c, c], [c, nc, c, c], [c, c, nc, c], [c, c, c, nc]]])
    def test_jc_q_mat(self):
        jc = JukesCantorModel()
        _LOG.debug("jc.q_mat = %s" % str(jc.q_mat))
        assert_mat_eq(self, jc.q_mat, [[-1.0, 1.0/3, 1.0/3, 1.0/3], [1.0/3, -1.0, 1.0/3, 1.0/3], [1.0/3, 1.0/3, -1.0, 1.0/3], [ 1.0/3, 1.0/3, 1.0/3, -1.0]])
    def test_rev_init(self):
        b = DNAType()
        a = RevDiscStateContTimeModel(r_upper=[[1.0, 1.0, 1.0], [1.0, 1.0], [1.0],], state_freq=(.25, .25, .25, .25, ), char_type=b)
        self.assertEqual(a.char_type, b)
    def test_initk2p(self):
        m = Kimura2ParameterModel(2.0)
        nc, ti, tv = 0.951679099289, 0.0239356129609, 0.0121926438748
        assert_list_of_mat_eq(self, m.calc_prob_matrices(0.05), [[[nc, tv, ti, tv], [tv, nc, tv, ti], [ti, tv, nc, tv], [tv, ti, tv, nc]]])
        self.assertRaises(TypeError, Kimura2ParameterModel)
        self.assertRaises(ValueError, Kimura2ParameterModel, 'hi there')
        
    def test_inithky(self):
        m = HKY85Model(2.0, [.25, 0.25, 0.25, 0.25])
        self.assertRaises(AttributeError, m.kappa.set_value, 3)
        nc, ti, tv = 0.951679099289, 0.0239356129609, 0.0121926438748
        assert_list_of_mat_eq(self, m.calc_prob_matrices(0.05), [[[nc, tv, ti, tv], [tv, nc, tv, ti], [ti, tv, nc, tv], [tv, ti, tv, nc]]])

    def test_bad(self):
        self.assertRaises(TypeError, RevDiscStateContTimeModel, [])
        self.assertRaises(ValueError, RevDiscStateContTimeModel, rmat=[[0,1],[.5,0]])

    def test_r_upper_to_r_mat(self):
        cases = [([[1.0]],
                  [[0.0, 1.0],
                   [1.0,0.0]]),
                 ([[1.0,1.0],
                   [1.0]],
                  [[0.0, 1.0, 1.0],
                   [1.0, 0.0, 1.0],
                   [1.0, 1.0, 0.0]]),
                 ([[1.0,2.0],
                   [3.0]],
                  [[0.0, 1.0, 2.0],
                   [1.0, 0.0, 3.0],
                   [2.0, 3.0, 0.0]]),
                ]
        for t_case in cases:
            rup, rmat = t_case
            r_mat_gen = _r_upper_to_r_mat(rup)
            r_up_gen = _r_mat_to_r_upper(rmat)
            assert_mat_eq(self, r_mat_gen, rmat)
            assert_mat_eq(self, r_up_gen, rup)
class Skip:
    def test_set_q_mat(self):
        a = RevDiscreteModel(r_upper=[[1.0, 1.0, 1.0], [1.0, 1.0], [1.0],])
        oth = 1.0/3.0
        assert_list_eq(self, a.state_freq, [0.25]*4)
        assert_mat_eq(self, a.q_mat, [[-1.0, oth, oth, oth],
                                       [oth, -1.0, oth, oth],
                                       [oth, oth, -1.0, oth],
                                       [oth, oth, oth, -1.0]])
        a.r_upper = [[1.0, 2.0, 1.0], [1.0, 2.0], [1.0],] # kimura kappa=2
        assert_list_eq(self, a.state_freq, [0.25]*4)
        assert_mat_eq(self, a.q_mat, [[-1.0, .25, .5, .25],
                                      [.25, -1.0, .25, .5],
                                      [.5, .25, -1.0, .25],
                                      [.25, .5, .25, -1.0]])
        a.state_freq = [0.2, 0.3, 0.15, 0.35]
        assert_list_eq(self, a.state_freq, [0.2, 0.3, 0.15, 0.35])
        exp = [[-0.9547738, 0.301507, 0.3015075, 0.351758],
               [0.20100502, -1.05527638, 0.15075376, 0.7035175],
               [0.40201005, 0.301507537, -1.05527638, 0.35175879],
               [0.2010050, 0.603015075, 0.150753768, -0.954773869]]
        assert_mat_eq(self, a.q_mat, exp)


def additional_tests():
    "returns all tests in this file as suite"
    return unittest.TestLoader().loadTestsFromTestCase(ModelTest)

# pylint: disable-msg=C0103
def getTestSuite():
    """Alias to the additional_tests().  This is unittest-style.
    `additional_tests` is used by setuptools.
    """
    return additional_tests()

if __name__ == "__main__":
    unittest.main()
    print "SKIPPING TESTS!!!"



##############################################################################
##  pytbeaglehon phylogenetic likelihood caluclations using beaglelib.
##
##  Copyright 2010 Mark T. Holder
##  All rights reserved.
##
##  See "LICENSE.txt" for terms and conditions of usage.
##
##############################################################################
