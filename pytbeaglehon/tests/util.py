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
    for ret_row, exp_row in izip(returned, expected):
        assert_list_eq(self, ret_row, exp_row)

def assert_list_eq(self, returned, expected):
    "calls self.assertAlmostEqual for all elements of two lists."
    for ret_val, exp_val in izip(returned, expected):
        self.assertAlmostEqual(ret_val, exp_val, places=5)

def additional_tests():
    "returns all tests in this file as suite"
    return unittest.TestSuite([])


def newick_tokenizer(newick):
    '''Takes a newick string serves as a generator for tokens.''' 
    import re
    r = re.compile(r'[(),:]')
    index = 0
    sp_tokens = r.split(newick)
    if (sp_tokens[0] != '') or (newick[0] != '('):
        raise ValueError("Expecting the string to start with '('")
    if (sp_tokens[-1] != '') or (newick[-1] != ')'):
        raise ValueError("Expecting the string to end with ')'")
    for t in sp_tokens[1:-1]:
        yield newick[index]
        if t.strip() != '':
            yield t.strip()
        index += 1 + len(t)
    yield newick[-1]
    
class NodeForTesting(object):
    def __init__(self, parent=None):
        self.parent = parent
        self.children = []
        self._edge_len = None
        self.leaf_index = None
    
    def parse(newick):
        curr_node = NodeForTesting()
        expect_edge_len =False
        for token in newick_tokenizer(newick):
            if expect_edge_len:
                curr_node.edge_len = float(token)
                expect_edge_len = False
            else:
                if token == '(':
                    n = NodeForTesting(parent=curr_node)
                    curr_node.children.append(n)
                    curr_node = n
                elif token == ')':
                    assert(curr_node.children or (curr_node.leaf_index is not None))
                    curr_node = curr_node.parent
                    assert(curr_node.leaf_index is None) # we don't handle internal taxon labels ...
                    if (curr_node is None):
                        raise ValueError("Newick error too many )")
                elif token == ',':
                    assert(curr_node.children or (curr_node.leaf_index is not None))
                    assert(curr_node.parent is not None)
                    n = NodeForTesting(parent=curr_node.parent)
                    curr_node.parent.children.append(n)
                    assert(curr_node.parent.leaf_index is None) # we don't handle internal taxon labels ...
                    curr_node = n
                elif token == ':':
                    expect_edge_len = True
                else:
                    try:
                        n = int(token) - 1
                        assert(n >= 0)
                    except:
                        raise ValueError("Expecting a (positive) taxon number but found %s" % token)
                    curr_node.leaf_index = n
        assert((curr_node.parent is None) and (expect_edge_len == False))
        return curr_node
    parse = staticmethod(parse)
class TreeForTesting(object):
    def __init__(self, newick):
        self.root = NodeForTesting.parse(newick)
        
    def postorder_internal_node_iter(self):
        st = []
        preorder = []
        curr_nd = self.root
        while True:
            c = curr_nd.children
            if len(c) == 0:
                try:
                    curr_nd = st.pop()
                except IndexError:
                    break
            else:
                preorder.append(curr_nd)
                assert(len(c) == 2)
                curr_nd = c[0]
                st.append(c[1])
            
        preorder.reverse()
        return preorder

##############################################################################
##  pytbeaglehon phylogenetic likelihood caluclations using beaglelib.
##
##  Copyright 2010 Mark T. Holder
##  All rights reserved.
##
##  See "LICENSE.txt" for terms and conditions of usage.
##
##############################################################################
