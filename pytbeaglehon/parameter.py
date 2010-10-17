#! /usr/bin/env python
'''Wrappers around numbers to encapsulate the concept of a free parameter of a 
model.
'''
# Copyright (c) 2010 by Mark T. Holder,  University of Kansas
# (see bottom of file)

class Parameter(object):
    def __init__(self, value):
        self._value = self.native_type()(value)
        self.listener_list = []
    def native_type(self):
        return lambda x : x
    def add_listener(self, listener):
        self.listener_list.append(listener)
    def get_value(self):
        return self._value
    def set_value(self, v):
        raise AttributeError("can't set attribute")
    value = property(get_value, set_value)

class FloatParameter(Parameter):
    def native_type(self):
        return float
    def __float__(self):
        return self.value
    def __add__(self, other):
        return float(self) + other
    def __sub__(self, other):
        return float(self) - other
    def __mul__(self, other):
        return float(self) * other
    def __div__(self, other):
        return float(self) / other
    def __repr__(self):
        return 'FloatParameter(%s)' % repr(self.value)

class MutableFloatParameter(FloatParameter):
    def set_value(self, v):
        self._value = v
    value = property(Parameter.get_value, set_value)
    

##############################################################################
##  pytbeaglehon phylogenetic likelihood caluclations using beaglelib.
##
##  Copyright 2010 Mark T. Holder
##  All rights reserved.
##
##  See "LICENSE.txt" for terms and conditions of usage.
##
##############################################################################
