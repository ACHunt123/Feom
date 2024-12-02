#!/usr/bin/env python
import numpy as np
import sys
import matplotlib.pyplot as plt

import heom.dvr as dvr
from heom.input import InputObj
from heom.type_dicts import potentials
from heom.general import hbar, pi, formatflt


file_input = sys.argv[1]
file_tcf_qq = "tcf_qq_exact"
inp = InputObj(file_input)
pot = potentials[inp.pot_type](inp)

q1 = inp.q1
q0 = inp.q0
dq = inp.dq
mass = inp.mass
beta = inp.beta
t0 = inp.t0
t1 = inp.t1
t_sample = inp.t_sample
omega = inp.tani_omega_a
gamma = inp.eta/inp.gamma
om = np.sqrt(omega**2 - (gamma**2)/4)

def make_classical_tcf(t):
    cf = (np.exp(-gamma*t/2)
            * (np.cos(om*t)+gamma*np.sin(om*t)/(2*om))
            / (beta*mass*omega**2))
    return cf

def matsubara(n):
    return 2*n*pi/(beta*hbar)

def gamma_term(t):
    n_lim = 20 # The number of Matsubara terms considered
    gm = 0
    for n in range(1,n_lim):
        omega_n = matsubara(n)
        gm += omega_n*np.exp(-omega_n*t)/(
                (omega**2+omega_n**2)**2 - (gamma**2)*(omega_n**2))
    gm *= (2*gamma/(2*mass*beta))
    return gm

def make_tcf(t):
    cf = -0*gamma_term(t) + (hbar*np.exp(-gamma*t/2)/(2*mass*om))*np.real(
            np.exp(1j*om*t)*np.tanh(beta*hbar*(om+1j*gamma/2)/2))
    #cf = (hbar*np.exp(-gamma*t/2)/(2*mass*om))*np.real(
    #        np.exp(1j*om*t)*np.tanh(beta*hbar*(om+1j*gamma/2)/2))
    return cf

ts = np.linspace(0,t1, int(t1/t_sample))
tcf = []
for i in range(len(ts)):
    tcf.append(make_tcf(ts[i]))

# Writing out
with open(file_tcf_qq, "w") as f:
    for i in range(len(ts)):
        f.write(formatflt(ts[i]))
        f.write(formatflt(tcf[i]))
        f.write("\n")

# Plotting
fig = plt.figure(figsize=(8, 6), dpi=160)
plt.plot(ts,tcf)
plt.title("QQ(t) TCF")
plt.show()
plt.clf()
plt.close()


