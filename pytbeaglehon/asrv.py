#! /usr/bin/env python
# Copyright (c) 2010 by Mark T. Holder
# (see bottom of file)
"""Provides RateHetManager class for calculating among-site rate variation."""
from pytbeaglehon.ccore.disc_state_cont_time_model import casrvo_ctor, casrvo_set_shape, casrvo_get_rates

_TOLERANCE = 1.0e-6
_MIN_GAMMA_SHAPE = 1.0e-07
_LOW_RATE_MIN_GAMMA_SHAPE = 1.0e-200

class RateHetType:
    "Enumeration of the kinds of among-site rate heterogeneity."
    # pylint: disable-msg=R0903,W0232
    GAMMA_EQ_CAT_MEDIAN = 0 # keep this 0 (or coordinate with the ccore.ASRVObj
    GAMMA_EQ_CAT_MEAN = 1
    LAST_GAMMA = 1
    ARBITRARY_RATES = 2
    MAX_RATE_HET_TYPE = 2

class RateHetManager(object):
    "Calculates the rate multipliers for different rate categories."
    def __init__(self,
                 rate_het_type=RateHetType.ARBITRARY_RATES,
                 **kwargs):
        """Currently just supports the discrete approximation of
        gamma-distributed rates.

        `rate_het_type` indicates the distribution used to calculate rates
            this should be a facet of the RateHetType "enum"

        """
        if rate_het_type > RateHetType.MAX_RATE_HET_TYPE:
            raise ValueError("Illegal value for rate_het_type")
        self._rate_het_type = rate_het_type
        num_categories = kwargs.get("num_categories", 0)
        rates = kwargs.get("rates")
        probabilities = kwargs.get("probabilities")
        if not num_categories:
            if rates:
                num_categories = len(rates)
            elif probabilities:
                num_categories = len(probabilities)
            else:
                raise ValueError('At least one of "num_categories", "rates", or "probabilities" arguments must be used')
        if num_categories < 1 or int(num_categories) != num_categories:
            raise ValueError("num_categories must be a positive integer")
        self._num_cat = num_categories
        self._probabilities = [1.0/num_categories] * num_categories
        self._rate_list = [1.0] * num_categories
        if rate_het_type <= RateHetType.LAST_GAMMA:
            shape = float(kwargs.get("shape", 0.5))
            if shape <= 0.0:
                raise ValueError("shape must be greater than 0.0")
            s = max(_MIN_GAMMA_SHAPE, shape)
            a = casrvo_ctor(s, num_categories, rate_het_type)
            casrvo_set_shape(a, shape)
            assert(a is not None)
            self._asrv = a
            self._rate_list = None
            self.shape = s

        if rates:
            self.rates = rates
        if probabilities:
            self.probabilities = probabilities


    def get_num_cat(self):
        "Returns the number of rate categories."
        return self._num_cat

    def get_rates(self):
        "Returns a list with a rate for each category."
        if not self._rate_list:
            if self._shape > _MIN_GAMMA_SHAPE:
                self._rate_list = casrvo_get_rates(self._asrv)
            else:
                self._rate_list = [_LOW_RATE_MIN_GAMMA_SHAPE ] * self._num_cat
                self._rate_list[self._num_cat - 1] = float(self._num_cat)
        return list(self._rate_list)

    def set_rates(self, r):
        "Returns a list with a rate for each category."
        if self._rate_het_type <= RateHetType.LAST_GAMMA:
            raise TypeError("RateHetManager.set_rates cannot be used with gamma distributions")
        if len(r) != self._num_cat:
            raise TypeError("Expecting a list of %d rates" % self._num_cat)
        rf = [float(i) for i in r]
        if min(rf) < 0.0:
            raise ValueError("All rates must be non-negative")
        self._rate_list = rf

    def get_probabilities(self):
        "Returns a list with a rate for each category."
        return list(self._probabilities)

    def set_probabilities(self, p):
        "Returns a list with a rate for each category."
        if self._rate_het_type <= RateHetType.LAST_GAMMA:
            raise TypeError("RateHetManager.set_probabilities cannot be used with gamma distributions")
        if len(p) != self._num_cat:
            raise TypeError("Expecting a list of %d probabilities" % self._num_cat)
        p = [float(i) for i in p]
        sp = sum(p)
        if sp < (1.0 - _TOLERANCE) or sp > (1.0 + _TOLERANCE):
            raise ValueError("Sum of probabilities must be 1.0")
        self._probabilities = p


    def get_shape(self):
        "Returns the shape parameter for a Gamma distribution over rates."
        if self._rate_het_type > RateHetType.LAST_GAMMA:
            raise TypeError("RateHetManager.get_shape can only be used with gamma distributions")
        return self._shape

    def set_shape(self, shape):
        "Sets the shape parameter for a Gamma distribution over rates."
        if self._rate_het_type > RateHetType.LAST_GAMMA:
            raise TypeError("RateHetManager.set_shape can only be used with gamma distributions")
        self._rate_list = None # this will trigger recalculation of the rates
        self._shape = shape
        if self._shape > _MIN_GAMMA_SHAPE:
            casrvo_set_shape(self._asrv, shape)

    num_categories = property(get_num_cat)
    rates = property(get_rates, set_rates)
    probabilities = property(get_probabilities, set_probabilities)
    shape = property(get_shape, set_shape)

def GammaRateHetManager(shape, num_categories, use_mean=True):
    """Factory for Gamma-distributed rates with equally frequent categories. If
        `use_mean` is False then the median rate will be used.
    """
    if use_mean:
        flag = RateHetType.GAMMA_EQ_CAT_MEAN
    else:
        flag = RateHetType.GAMMA_EQ_CAT_MEDIAN
    return RateHetManager(shape=shape, num_categories=num_categories, rate_het_type=flag)

##############################################################################
##  pytbeaglehon phylogenetic likelihood caluclations using beaglelib.
##
##  Copyright 2010 Mark T. Holder
##  All rights reserved.
##
##  See "LICENSE.txt" for terms and conditions of usage.
##
##############################################################################
