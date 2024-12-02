#!/usr/bin/env python
# File: dvr.py
"""DVR derivative matrices"""

import numpy as np
from heom.general import hbar

def derivative1(xs):
    """Makes a first derivative matrix."""
    n = len(xs)+1
    L = xs[-1]-xs[0]
    D = np.zeros([n-1,n-1])

    def fun_sum(A,n):
        return -((1-n)*np.sin(n*A)+n*np.sin((n-1)*A))/(2*(1-np.cos(A)))

    for i in range(1,n):
        for j in range(1,n):
            A = (i-j)*np.pi/n
            B = (i+j)*np.pi/n
            if i==j:
                D[i-1,j-1] = (np.pi/(n*L))*fun_sum(B,n)
            else:
                D[i-1,j-1] = (np.pi/(n*L))*(fun_sum(A,n)+fun_sum(B,n))
    return D

def derivative2(xs, dx):
    """Makes a second derivative matrix."""
    N = len(xs)
    T = np.zeros([N,N])
    for i in range(1,N+1):
        for j in range(1, N+1):
            if i==j:
                T[i-1,j-1] = -( (-1)**(i-j) /(dx**2) )*np.pi**2 /3
            else:
                T[i-1,j-1] = -( (-1)**(i-j) /(dx**2) )*2/((i-j)**2)
    return T

def derivative2_finite(xs, dx):
    """Makes a second derivative matrix."""
    N = len(xs)+1
    n = N-1
    T = np.zeros([n,n])
    a = xs[0]-dx
    b = xs[-1]+dx
    for i in range(1,N):
        for j in range(1, N):
            if i==j:
                T[i-1,j-1] = -( 1 /((b-a)**2) ) * (np.pi**2 /2) * ( (2*N*N+1)/3 - 1/(np.sin(np.pi*i/N)**2))
            else:
                T[i-1,j-1] = -( (-1)**(i-j) /((b-a)**2) ) * (np.pi**2 /2) * ( 1/(np.sin(np.pi*(i-j)/(2*N))**2) - 1/(np.sin(np.pi*(i+j)/(2*N))**2))
    return T

def hamiltonian(qs, pot, mass):
    n_q = len(qs)
    dq = (qs[-1]-qs[0])/(n_q-1)
    kinE = -(hbar**2/(2*mass))*derivative2(qs,dq)
    potE = np.identity(n_q)
    for i in range(n_q):
        potE[i,i] = pot.calc(qs[i])
    ham = kinE + potE
    return ham

def hamiltonian_finite(qs, pot, mass):
    n_q = len(qs)
    dq = (qs[-1]-qs[0])/(n_q-1)
    kinE = -(hbar**2/(2*mass))*derivative2_finite(qs,dq)
    potE = np.identity(n_q)
    for i in range(n_q):
        potE[i,i] = pot.calc(qs[i])
    ham = kinE + potE
    return ham

def ham_evals_evecs(qs, pot, mass):
    n_q = len(qs)
    dq = (qs[-1]-qs[0])/(n_q-1)
    ham = hamiltonian(qs, pot, mass)
    evals, evecs = np.linalg.eigh(ham)
    # The function returns the transformation matrix of normalised eigenvectors
    #|  |    |    |  |
    #|  |    |    |  |
    #|  |    |    |  |
    #|psi1 psi2 psi3 |
    #|  |    |    |  |
    #|  |    |    |  |
    #|  |    |    |  |
    # These wavefunction are in fact psi(q_i)*sqrt(dq), not the actual values of psi(q_i)
    # This means that the normalisation condition is psi dot psi = 1
    return ham, evals, evecs
