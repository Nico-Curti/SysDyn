#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sympy as sy

__package__ = "Lotka Volterra"
__author__  = "Nico Curti"
__email__   = "nico.curti2@unibo.it"

def lotka_volterra(y, t, A, B, C, D):
  x1, x2 = y
  return (A - B * x2) * x1, (C * x1 - D) * x2


if __name__ == '__main__':

  A, B, C, D = sy.symbols(['A', 'B', 'C', 'D'], positive=True)
  x, y = sy.symbols(['x', 'y'], positive=True)

  dx = sy.Eq(x, (A - B * y) * x )
  dy = sy.Eq(y, (C * x - D) * y )

  solutions = sy.solve([dx, dy], [x, y])

  print(solutions)
  print(sy.latex(solutions))

