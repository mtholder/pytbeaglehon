#! /usr/bin/env python
# Copyright (c) 2005-7 by Mark T. Holder,  University of Kansas
# (see bottom of file)

"unit tests datatyp internals"
import unittest
from pytbeaglehon.tests.util import *
# pylint: disable-msg=C0111,W0401,W0611,W0212
from pytbeaglehon.disc_char_type import DiscreteCharType, DNAType, AAType, AANoStopType

class DiscreteCharTypeTest(unittest.TestCase):
    def test_DNA(self):
        d = DNAType()
        inds = d.to_indices("ACNGT-WAYKCSXBDVMRH")
        self.assertEquals(inds, [0, 1, 4, 2, 3, 4, 8, 0, 6, 10, 1, 9, 4, 14, 13, 11, 7, 5, 12])
        self.assertEquals(d.num_states, 4)
        self.assertEquals(d.states, ('A', 'C', 'G', 'T'))
        self.assertEquals(d.all_symbols, ('A', 'C', 'G', 'T', 'N', 'R', 'Y', 'M', 'W', 'S', 'K', 'V', 'H', 'D', 'B'))
        self.assertEquals(d.ignore_case, True)
        self.assertEquals(d.symbol_to_ind, {'A': 0, 'C': 1, 'B': 14, 'D': 13, 'G': 2, 'H': 12, 'K': 10, 'M': 7, 'N': 4, 'S': 9, 'R': 5, 'T': 3, 'W': 8, 'V': 11, 'Y': 6, 'X': 4, '-': 4, '?': 4})
        self.assertEquals(d.state_code_lookup, ((0, ), (1, ), (2, ), (3, ), (0, 1, 2, 3), (0, 2), (1, 3), (0, 1), (0, 3), (1, 2), (2, 3), (0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)))
        self.assertTrue(d.state_code_lookup is not None)
        self.assertEquals(d.partial_ambiguity_indices, tuple(range(5,15)))
        self.assertFalse(d.has_partial_ambiguity([1, 2, 4 , 0, 0, 0, 0, 2, 3, 4]))
        for i in range(5, 15):
            self.assertTrue(d.has_partial_ambiguity([1, 2, 4, 0, i, 0, 3, 0, 2, 3, 4]))

    def test_AA(self):
        d = AAType()
        inds = d.to_indices("ACDEFGHIKLMNPQRSTVWY*XBZ?-")
        self.assertEquals(inds, range(24) + [21, 21])
        self.assertEquals(d.num_states, 21)
        self.assertEquals(d.states, ('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', '*'))
        self.assertEquals(d.all_symbols, ('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', '*', 'X', 'B', 'Z'))
        self.assertEquals(d.ignore_case, True)
        self.assertEquals(d.symbol_to_ind, {'A': 0, 'C': 1, 'B': 22, 'E': 3, 'D': 2, 'G': 5, 'F': 4, 'I': 7, 'H': 6, 'K': 8, 'M': 10, 'L': 9, 'N': 11, 'Q': 13, 'P': 12, 'S': 15, 'R': 14, 'T': 16, 'W': 18, 'V': 17, 'Y': 19, 'X': 21, 'Z': 23, '*': 20, '-': 21, '?': 21})
        self.assertEquals(d.state_code_lookup, ((0,), (1,), (2,), (3,), (4,), (5,), (6,), (7,), (8,), (9,), (10,), (11,), (12,), (13,), (14,), (15,), (16,), (17,), (18,), (19,), (20,), (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20), (2, 11), (3, 13)))
        self.assertTrue(d.state_code_lookup is not None)
        self.assertEquals(d.partial_ambiguity_indices, tuple(range(22,24)))
        self.assertFalse(d.has_partial_ambiguity([1, 2, 4 , 0, 0, 0, 0, 2, 3, 4]))
        for i in [22, 23]:
            self.assertTrue(d.has_partial_ambiguity([1, 2, 4, 0, i, 0, 3, 0, 2, 3, 4] + range(21)))

    def test_AANoStop(self):
        d = AANoStopType()
        inds = d.to_indices("ACDEFGHIKLMNPQRSTVWYXBZ?-")
        self.assertEquals(inds, range(23) + [20, 20])
        self.assertEquals(d.num_states, 20)
        self.assertEquals(d.states, ('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'))
        self.assertEquals(d.all_symbols, ('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', 'X', 'B', 'Z'))
        self.assertEquals(d.ignore_case, True)
        self.assertEquals(d.symbol_to_ind, {'A': 0, 'C': 1, 'B': 21, 'E': 3, 'D': 2, 'G': 5, 'F': 4, 'I': 7, 'H': 6, 'K': 8, 'M': 10, 'L': 9, 'N': 11, 'Q': 13, 'P': 12, 'S': 15, 'R': 14, 'T': 16, 'W': 18, 'V': 17, 'Y': 19, 'X': 20, 'Z': 22, '-': 20, '?': 20})
        self.assertEquals(d.state_code_lookup, ((0,), (1,), (2,), (3,), (4,), (5,), (6,), (7,), (8,), (9,), (10,), (11,), (12,), (13,), (14,), (15,), (16,), (17,), (18,), (19,), (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19), (2, 11), (3, 13)))
        self.assertTrue(d.state_code_lookup is not None)
        self.assertEquals(d.partial_ambiguity_indices, tuple(range(21, 23)))
        self.assertFalse(d.has_partial_ambiguity([1, 2, 4 , 0, 0, 0, 0, 2, 3, 4]))
        for i in [21, 22]:
            self.assertTrue(d.has_partial_ambiguity([1, 2, 4, 0, i, 0, 3, 0, 2, 3, 4] + range(20)))


    def test_bad(self):
        self.assertRaises(ValueError, DiscreteCharType, "")
        self.assertRaises(ValueError, DiscreteCharType, "ACA")
        self.assertRaises(ValueError, DiscreteCharType, "ACGT", (("A", "A"), ))
        self.assertRaises(ValueError, DiscreteCharType, "ACGT", (("W", "AT"), ("W", "AT")))
        self.assertRaises(ValueError, DiscreteCharType, "ACGT", (("W", "AT"), ("Y", "CK")))
        self.assertRaises(ValueError, DiscreteCharType, "ACGT", (("W", "AT"), ("Y", "CT")), missing="W")
        self.assertRaises(ValueError, DiscreteCharType, "ACGT", (("W", "AT"), ("Y", "CT")), missing="?", aliases=[("A", "A")])
        self.assertRaises(ValueError, DiscreteCharType, "ACGT", (("W", "AT"), ("Y", "CT")), missing="?", aliases=[("a", "a")], ignore_case=False)

def additional_tests():
    "returns all tests in this file as suite"
    return unittest.TestLoader().loadTestsFromTestCase(DiscreteCharTypeTest)

# pylint: disable-msg=C0103
def getTestSuite():
    """Alias to the additional_tests().  This is unittest-style.
    `additional_tests` is used by setuptools.
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
