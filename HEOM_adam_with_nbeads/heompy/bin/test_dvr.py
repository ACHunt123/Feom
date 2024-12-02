#!/usr/bin/env python
# File: test_dvr.py
import numpy as np

import heom.dvr as dvr
import matplotlib.pyplot as plt

span = 10
npoints = 1000
xmin = -span
xmax = span
xs = np.linspace(xmin,xmax,npoints)
dx = xs[1]-xs[0]

def f(x):
    return (-x**2 + x**4)*np.exp(-x**2)
def df(x):
    return (-2*x + 4*x**3)*np.exp(-x**2) + (-x**2 + x**4)*np.exp(-x**2)*(-2*x)
def df2(x):
    return np.exp(-x**2)*(4*x**6 - 22*x**4 + 22 * x**2 -2)

ys0 = f(xs)
# Analytic derivative
ys1 = df(xs)
# Numerical derivative
mat = dvr.derivative1(xs)
ys1n = mat @ ys0

ys2 = df2(xs)
mat2 = dvr.derivative2(xs,dx)
ys2n = mat2 @ ys0

plt.plot(xs, ys1n, label="f' numerical")
plt.plot(xs, ys1, label="f' analytic", ls='--')
plt.plot(xs, ys2n, label="f'' numerical")
plt.plot(xs, ys2, label="f'' analytic", ls='--')
plt.legend()
plt.show()
