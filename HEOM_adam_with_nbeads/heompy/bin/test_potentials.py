#!/usr/bin/env python
# File: pytest_potentials.py
import os
import numpy as np
import matplotlib.pyplot as plt

from heom.input import InputObj
from heom.type_dicts import potentials

filename = os.environ['WORKDIR']+ "/heompy/input/input"
inp = InputObj(filename)
print(inp.__dict__)

xs = np.linspace(-10, 10, 1000)

plt.ylim([-100, 500])
for i in range(len(potentials)):
    pot = potentials[i](inp)
    ys = pot.calc(xs)
    plt.plot(xs, ys, label=pot.__class__.__name__)
plt.title("Potentials")
plt.legend()
plt.show()
plt.clf()

plt.ylim([-500, 500])
for i in range(len(potentials)):
    pot = potentials[i](inp)
    ys = pot.diff(xs)
    plt.plot(xs, ys, label=pot.__class__.__name__)
plt.title("Derivatives")
plt.legend()
plt.show()

plt.ylim([-100, 500])
for i in range(len(potentials)):
    pot = potentials[i](inp)
    pot.switch()
    ys = pot.calc(xs)
    plt.plot(xs, ys, label=pot.__class__.__name__)
plt.title("Potentials - Switch")
plt.legend()
plt.show()
plt.clf()

plt.ylim([-500, 500])
for i in range(len(potentials)):
    pot = potentials[i](inp)
    pot.switch()
    ys = pot.diff(xs)
    plt.plot(xs, ys, label=pot.__class__.__name__)
plt.title("Derivatives - Switch")
plt.legend()
plt.show()
