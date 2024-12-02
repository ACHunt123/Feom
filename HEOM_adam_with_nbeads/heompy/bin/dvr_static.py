#!/usr/bin/env python
import numpy as np
import sys
import matplotlib.pyplot as plt

import heom.dvr as dvr
#from heom.input import InputObj
#from heom.type_dicts import potentials, states
from heom.general import hbar, pi, formatflt

class PotBox(object):
    """Harmonic oscillator potential"""
    def __init__(self):
        pass
    def calc(self,x):
        return 0

class PotHO(object):
    """Harmonic oscillator potential"""
    def __init__(self):
        self.omega = 1
        self.mass = 1
        self.k = self.mass*(self.omega**2)

    def calc(self,x):
        return 0.5*self.k*x**2

pot_dict = {0: PotBox, 1: PotHO}

pot_type = 1
pot = pot_dict[pot_type]()
mass = 1.0

rang = 3
a = -rang
b = rang
dx = 0.1
xs = np.arange(a,b,dx)
ham = dvr.hamiltonian_finite(xs, pot, mass)
ham2 = dvr.hamiltonian(xs, pot, mass)

evals, evecs = np.linalg.eigh(ham)
evals2, evecs2 = np.linalg.eigh(ham2)
evalnum = len(evals)
evalnum2 = len(evals2)
# Normalising Wavefunctions
for n in range(evalnum):
    evecs[:,n] /= np.sqrt(dx*np.vdot(evecs[:,n],evecs[:,n])) 
for n in range(evalnum2):
    evecs2[:,n] /= np.sqrt(dx*np.vdot(evecs2[:,n],evecs2[:,n])) 


print("{:>19s}{:>19s}{:>19s}{:>19s}".format("number","exact","computed finite","computed"))
for ind in range(evalnum):
    if pot_type == 0:
        QN = ind+1
        L = b-a+2*dx
        E_exact = np.pi**2 * QN**2 / (2*mass*L**2)
    elif pot_type == 1:
        QN = ind
        E_exact = (QN + 0.5)
    print("{:19d}{:19f}{:19f}{:19f}".format(ind,E_exact,evals[ind],evals2[ind]))
