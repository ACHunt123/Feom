#!/usr/bin/env python3
# File: psd.py
import numpy as np

scheme = "pade_N/N"
scheme = "pade_N+1/N"
scheme = "pade_N-1/N"

def sgn(x):
    if x>0:
        return 1
    elif x<0:
        return -1
    else:
        return 0


def delta(i,j):
    if i==j:
        return 1
    else:
        return 0

def b(n):
    return 2*n+1

def d(n):
    if n==1:
        return 1/b(1)
    elif n%2==0:
        m = n//2
        return -4 * m**2 * b(m)**2 * b(2*m)
    elif n%2==1:
        m = (n-1)//2
        return -b(2*m+1)/(4*m*(m+1)*b(m)*b(m+1))

def Tv(N):
    s = 0
    for n in range(1,N+2):
        s += d(2*n)
    return 1/(4*s)

def Rv(N):
    s = 0
    for n in range(1,N+2):
        s += d(2*n-1)
    s2 = 0
    for n in range(1,N+2):
        s += d(2*n)
    return (4*Tv(N))**2 * s**2 * s2**2

def t(N):
    if N==1:
        return Tv(0)
    else:
        return Tv(N-1)/Tv(N-2)

def get_Lambda(N):
    if scheme=="pade_N-1/N":
        limit = 2*N
    elif scheme in {"pade_N/N","pade_N+1/N"}:
        limit = 2*N+1
    Lambda = np.zeros([limit,limit])
    for i1 in range(limit):
        for i2 in range(limit):
            m = i1+1
            n = i2+1
            if scheme in {"pade_N-1/N","pade_N/N"}:
                Lambda[i1,i2] = (delta(m,n+1)+delta(m,n-1))/np.sqrt(b(m)*b(n))
            elif scheme=="pade_N+1/N":
                Lambda[i1,i2] = (delta(m,n+1)+delta(m,n-1))/np.sqrt(d(m+1)*d(n+1))
    return Lambda

def get_Lambda_tilde(N):
    if scheme in {"pade_N-1/N","pade_N/N"}:
        Lambda = get_Lambda(N)[1:,1:]
    elif scheme=="pade_N+1/N":
        limit = 2*N+2
        Lambda = np.zeros([limit,limit],dtype=complex)
        for i1 in range(limit):
            for i2 in range(limit):
                m = i1+1
                n = i2+1
                Lambda[i1,i2] = (delta(m,n+1)+delta(m,n-1))/np.sqrt(d(m)*d(n))
    return Lambda

def get_xi_zeta(N):
    evals = np.linalg.eigvalsh(get_Lambda(N))
    xi = []
    for e in evals:
        if e > 1e-16:
            xi.append(2/e)
    zeta = []
    if scheme in {"pade_N-1/N","pade_N/N"}:
        evals = np.linalg.eigvalsh(get_Lambda_tilde(N)) # i.e. Lambda tilde
        for e in evals:
            if e > 1e-16:
                zeta.append(2/e)
    return np.sort(np.array(xi)),np.sort(np.array(zeta))

def get_eta(N,zeta,xi):
    eta = np.zeros([N])
    if scheme=="pade_N-1/N":
        for j in range(N):
            eta[j] = (1/2)*N*b(N+1)
            for k in range(N):
                if k<(N-1):
                    eta[j] *= zeta[k]**2 - xi[j]**2
                if k!=j:
                    eta[j] /= xi[k]**2 - xi[j]**2
    elif scheme=="pade_N/N":
        for j in range(N):
            Rn = 1/(4*(N+1)*b(N+1))
            eta[j] = (1/2)*Rn
            for k in range(N):
                eta[j] *= zeta[k]**2 - xi[j]**2
                if k!=j:
                    eta[j] /= xi[k]**2 - xi[j]**2
    elif scheme=="pade_N+1/N":
        #*************
        def i2m(i):
            return i-1
        def m2i(m):
            return m+1
        #*************
        for j_ind in range(N):
            j=j_ind+1
            # Making deltas and rs
            def delta_diff(k):
                k_ind = k-1
                # note that k is the index k = real k-1
                if k==j or k==N+1:
                    return 1
                else:
                    return xi[k_ind]**2-xi[j_ind]**2
            def r(kin):
                if kin==0:
                    return 0
                elif kin%2==1:
                    k = (kin+1)//2
                    frac = t(k)/delta_diff(k)
                    return np.sqrt(4*abs(frac))
                else:
                    k = kin//2
                    frac = t(k)/delta_diff(k)
                    return np.sqrt(4*abs(frac))*sgn(frac)
            # *******************
            Xs = np.zeros([2*N+4])
            Xs[m2i(-1)] = 0
            Xs[m2i(0)]  = 1/2
            for m in range(1,2*N+3):
                Xs[m2i(m)] = r(m)*(d(m)*Xs[m2i(m-1)]-r(m-1)* xi[j_ind]**2 * Xs[m2i(m-2)]/4)
            eta[j_ind] = Xs[m2i(2*N+2)]
    return eta

def get_xi_eta(N):
    xi,zeta = get_xi_zeta(N)
    eta = get_eta(N,zeta,xi)
    return xi,eta

def get_exact(x):
    return 1/(1-np.exp(-x))

def get_PSD(N):
    xi,eta = get_xi_eta(N)
    if scheme=="pade_N-1/N":
        def fun(x):
            Phi = 0
            for j in range(N):
                Phi += 2*eta[j]/(x**2+xi[j]**2)
            return 1/2 + 1/x + x*Phi
    elif scheme=="pade_N/N":
        def fun(x):
            Phi = x/(4*(N+1)*b(N+1)) # x*Rn
            for j in range(N):
                Phi += 2*eta[j]/(x**2+xi[j]**2)
            return 1/2 + 1/x + x*Phi
    elif scheme=="pade_N+1/N":
        def fun(x):
            Phi = Rv(N) + Tv(N)*x**2 # x*Rn
            for j in range(N):
                Phi += 2*eta[j]/(x**2+xi[j]**2)
            return 1/2 + 1/x + x*Phi
    return fun

