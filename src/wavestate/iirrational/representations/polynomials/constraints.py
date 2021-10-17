# -*- coding: utf-8 -*-
"""
"""



class PolyConstraints(object):
    """
    Constraint flags used for polynomials. These have correspondences
    for the root constraints for most polynomial classes
    """
    no_constraint = frozenset()
    even_real   = frozenset(['even_real'])
    odd_real    = frozenset(['odd_real'])
    odd_imag    = frozenset(['odd_imag'])
    palendromic = frozenset(['palendromic'])
    odd_zero    = odd_imag | odd_real
    eRoR        = even_real | odd_real
    eRoI        = even_real | odd_imag
    eRoZ        = even_real | odd_zero
    eRoRpal     = even_real | odd_real | palendromic
    eRoIpal     = even_real | odd_imag | palendromic
    eRoZpal     = even_real | odd_zero | palendromic

poly_constraints = PolyConstraints()

