#! /usr/bin/env python
'''Wrappers around numbers to encapsulate the concept of a free parameter of a 
model.
'''
# Copyright (c) 2010 by Mark T. Holder,  University of Kansas
# (see bottom of file)

class Parameter(object):
    def __init__(self, value, **kwargs):
        self._prev_value = None
        self._value = self.native_type()(value)
        self.listener_list = []
        self.name = kwargs.get('name')
        self.is_mutable = False
        self._min_value = kwargs.get('min_value')
        self._max_value = kwargs.get('max_value')
        self._parent_param = kwargs.get('parent_param')
    def get_min_value(self):
        return self._min_value
    def set_min_value(self, x):
        if self.value < x:
            self._value = x
        self._min_value = x
    min_value = property(get_min_value, set_min_value)

    def get_max_value(self):
        return self._max_value
    def set_max_value(self, x):
        if self.value > x:
            self._value = x
        self._max_value = x
    max_value = property(get_max_value, set_max_value)
    def native_type(self):
        return lambda x : x
    def add_listener(self, listener):
        self.listener_list.append(listener)
    def del_listeneter(self, listener):
        try:
            self.listener_list.remove(listener)
            return True
        except:
            return False
    def get_value(self):
        return self._value
    def set_value(self, v, notify=True, notify_parent=True):
        if self.is_mutable:
            pv = self._prev_value
            if (self._min_value is not None) and (v < self._min_value):
                raise ValueError("Value of parameter [%s] must be >= %s" % (str(self), str(self._min_value)))
            if (self._max_value is not None) and (v > self._max_value):
                raise ValueError("Value of parameter [%s] must be <= %s" % (str(self), str(self._max_value)))
            to_set = self.native_type()(v)
            self._prev_value = self._value
            self._value = to_set
            if notify:
                try:
                    for listener in self.listener_list:
                        listener(self)
                except Exception, x:
                    self._value = _prev_value
                    self._prev_value = pv
                    try:
                        for listener in self.listener_list:
                            listener(self)
                    except:
                        pass
                    raise x
            if notify_parent and self._parent_param:
                self._parent_param.sub_parameter_changed(self)
        else:
            raise AttributeError("can't set attribute")
    value = property(get_value, set_value)
    def __str__(self):
        i = self.name and (' named "%s"' % self.name) or ''
        return 'Parameter(%s) at %d' % (str(self._value),i, id(self))
    def parameters(self):
        return [self]

class FloatParameter(Parameter):
    def __init__(self, value, **kwargs):
        Parameter.__init__(self, value, **kwargs)
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
    def __str__(self):
        i = self.name and (' named "%s"' % self.name) or ''
        return 'FloatParameter(%f)%s at %d' % (self._value, i,  id(self))

class MutableFloatParameter(FloatParameter):
    def __init__(self, value, **kwargs):
        FloatParameter.__init__(self, value, **kwargs)
        self.is_mutable = True
    def __str__(self):
        i = self.name and (' named "%s"' % self.name) or ''
        return 'MutableFloatParameter(%f)%s at %d' % (self._value, i, id(self))

class ProbabilityVectorParameter(Parameter):
    def __init__(self, elements, **kwargs):
        self._vec_len = len(elements)
        if self._vec_len < 2:
            raise ValueError("A ProbabilityVectorParameter must be off length 2 or greater")
        sub_p = []
        self._fixed_sub_p = set()
        self._free_parameters = set()
        all_p = True
        s = 0.0
        for el in elements:
            if isinstance(el, Parameter):
                if not el.is_mutable:
                    self._fixed_sub_p
                    self._fixed_sub_p.add(el)
                else:
                    self._free_parameters.add(el)
                sub_p.append(el)
                s += el.value
            else:
                np = FloatParameter(el)
                self._fixed_sub_p.add(np)
                sub_p.append(np)
                s += np.value
        if abs(1.0 - s > 1.0e-6):
            raise ValueError("Values in ProbabilityVectorParameter are expected to add to 1.0")
        if all_p and (len(self._fixed_sub_p) == (self._vec_len - 1)) :
            raise ValueError("Cannot have all but one parameter be fixed")
        if len(self._fixed_sub_p) < (self._vec_len - 1):
            self.is_mutable = True
        for n, el in enumerate(sub_p):
            el.max_value = 1.0
            el.min_value = 0.0
            el._parent_param = self
            el.index = n
        self.sub_parameters = sub_p

    def sub_parameter_changed(self, s):
        assert(s in self.sub_parameters)
        others = [i for i in self._free_parameters if i is not s]
        assert(len(others) > 0)
        other_prev_sum = sum([i.value for i in others])
        s_diff = s.value - s._prev_value
        if s_diff > other_prev_sum:
            raise ValueError("Setting parameter [%s] to %f forces the probability vector to exceed 1.0" % (str(s), s._value))
        if s_diff != 0.0:
            o_factor = (other_prev_sum - s_diff)/other_prev_sum
            for o in others:
                n = o._value * o_factor
                o.set_value(n, notify_parent=False)
            

        

##############################################################################
##  pytbeaglehon phylogenetic likelihood caluclations using beaglelib.
##
##  Copyright 2010 Mark T. Holder
##  All rights reserved.
##
##  See "LICENSE.txt" for terms and conditions of usage.
##
##############################################################################
