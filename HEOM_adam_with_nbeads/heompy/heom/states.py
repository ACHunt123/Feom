#!/usr/bin/env python
# File: states.py
""" Initial states for simulations."""
import numpy as np
import pickle
import os
import copy
from math import sinh, cosh, sqrt

from heom.general import pi, hbar
import heom.dvr as dvr
#from mpl_toolkits.mplot3d import Axes3D
#import matplotlib as mpl
#import matplotlib.pyplot as plt

#******************************************************************************
#******************************************************************************
def rho2wig(rho,qs):
    """Carries out numerical Wigner transform

    W(q,p) = const int <q-y/2|rho|q+y/2> exp(ipy/hbar) dy
    W(q,p) = const int <q-z|rho|q+z> exp(ip2z/hbar) dz
    W(q,p) = const int (<q-z|rho|q+z>)^* exp(-ip2z/hbar) dz
    where we have changed the variables and then used the fact
    that Wigner function W is real (W^* = W)

    Fourier transform as defined in NumPy is:
    FT(y) = const int f(x) exp(-i 2pi x y)

    The function first calculates the matrix elements <q-z|rho|q+z>
    by extending the density matrix so that the range of z can be kept
    constant. Then it carries out the Fourier transform, carefully
    using (i)fftshift, to get the correct answer.

    For the frequencies:
    - normal FT: I(nu) = ... exp(-i z 2 pi nu)
    - my FT:     I(p)  = ... exp(-i z 2p/hbar)
    Thus p = pi hbar nu,
    where nu is the frequency obtained from FT
    """
    n_q = len(qs)
    dq = (qs[-1]-qs[0])/(n_q-1)
    span = qs[-1]-qs[0]
    # Extending rho and qs
    rho_ext = np.zeros([3*n_q,3*n_q], dtype=complex)
    rho_ext[n_q:2*n_q,n_q:2*n_q] = rho
    qs_ext = np.zeros([3*n_q])
    qs_ext[n_q:2*n_q] = qs
    for i in range(n_q):
        qs_ext[i] = qs[i] - (span+dq)
        qs_ext[2*n_q+i] = qs[i] + (span+dq)
    # Calculating matrix elements f(q,z) = <q-z|rho|q+z> 
    fs = np.zeros([n_q,n_q], dtype=complex) # rows = qs, columns = zs
    iz_min = (-(n_q-1))//2
    iz_max = (n_q+1)//2 # iz reaches only iz_max-1
    for iq in range(n_q,2*n_q):
        for iz in range(iz_min,iz_max):
            # for odd n_q, iz goes from -(N-1)/2 to (N-1)/2
            # for even n_q, iz goes from -N/2 to (N-1)/2
            ket = np.zeros([3*n_q])
            ket[iq+iz] = 1
            bra = np.zeros([3*n_q])
            bra[iq-iz] = 1
            fs[iq-n_q,iz-iz_min] = bra@rho_ext@ket

    fs = np.conj(fs)

    # Numerical Wigner transform
    for iq in range(n_q):
        # Changing z interval from -T/2, T/2 to 0, T
        fs[iq,:] = np.fft.ifftshift(fs[iq,:])
    wig = np.zeros([n_q,n_q])
    for iq in range(n_q):
        # Fourier transform + appropriate shift in frequencies
        wig[iq,:] = np.fft.fftshift(np.fft.hfft(fs[iq,:],n=n_q))
    freqs =  np.fft.fftshift(np.fft.fftfreq(n_q,dq))
    ps = np.pi*hbar*freqs
    wig = np.real(wig)
    return (wig, ps)

#******************************************************************************
#******************************************************************************

class StateBase(object):
    """Initial state base class"""
    def __init__(self,inp,pot):
        self.mass = inp.mass
        self.pot = copy.deepcopy(pot)
        self.beta = inp.beta

    def rho(self, qs):
        raise NotImplementedError(
                "rho not implemented for this initial state.")

    def wigner(self, qs, dq):
        # Using numerical Wigner transform
        wig, ps = rho2wig(self.rho(qs),qs)
        dp = (ps[-1]-ps[0])/(len(ps)-1)

        return wig, ps, dp

    def rhoq2rhoe(self,rho,qs):
        #fig = plt.figure()
        #ax = fig.gca(projection='3d')
        #X,Y = np.meshgrid(qs,qs)
        #surf = ax.plot_surface(X,Y, rho.real, antialiased=True, cmap=mpl.cm.coolwarm)
        #plt.show()
        #plt.clf()
        print("Converting rhoq to rhoe using a potential that has a_ren",
                self.pot.a_ren)
        ham, evals, evecs = dvr.ham_evals_evecs(qs, self.pot, self.mass)
        dq = (qs[-1]-qs[0])/(len(qs)-1)
        rhoe = evecs.T.conj()@rho@evecs
        #fig = plt.figure()
        #ax = fig.gca(projection='3d')
        #X,Y = np.meshgrid(qs,qs)
        #surf = ax.plot_surface(X,Y, rhoe.real, antialiased=True, cmap=mpl.cm.coolwarm)
        #plt.show()
        #plt.clf()
        return rhoe, evals, evecs

    def rho_EE(self, qs, n_EE):
        """Generates energy eigenvector representation

        Can always be done by taking the position representation and
        transforming it to the energy representation.
        """
        rho, evals, evecs = self.rhoq2rhoe(self.rho(qs),qs)
        return rho[:n_EE,:n_EE], evals[:], evecs[:,:]
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

    def rho_q_kubo_q(self,qs):
        raise NotImplementedError(
                "rho_q_kubo_q not implemented for this initial state.")

    def rho_EE_kubo_q(self, qs, n_EE, q_mat):
        rho, evals, evecs = self.rhoq2rhoe(self.rho_q_kubo_q(qs),qs)
        return rho[:n_EE,:n_EE], evals[:], evecs[:,:]

    def rho_q_kubo_q2(self,qs):
        raise NotImplementedError(
                "rho_q_kubo_q not implemented for this initial state.")

    def rho_EE_kubo_q2(self, qs, n_EE, q_mat):
        rho, evals, evecs = self.rhoq2rhoe(self.rho_q_kubo_q2(qs),qs)
        return rho[:n_EE,:n_EE], evals[:], evecs[:,:]
#******************************************************************************
#******************************************************************************

class StateWP(StateBase):
    """Gaussian wave packet initial state."""
    def __init__(self, inp, pot):
        StateBase.__init__(self,inp,pot)
        self.q0 = inp.wp_q0
        self.p0 = inp.wp_p0
        self.sigma = inp.wp_sigma
        self.var = self.sigma**2

    def rho(self, qs):
        n_q = len(qs)
        dq = (qs[-1]-qs[0])/(n_q-1)
        psi = ((2*pi*self.var)**(-0.25)
                * np.exp(-((qs-self.q0)**2)/(4*self.var))
                * np.exp(1j*self.p0*qs/hbar) )
        rho0 = np.empty([n_q,n_q], dtype=complex)
        rho0[:,:] = np.outer(psi,psi)*dq
        return rho0

#******************************************************************************
#******************************************************************************

class StateTani(StateBase):
    """Harmonic oscillator in thermal equilibrium initial state. (Tanimura init. state)"""
    def __init__(self,inp,pot):
        StateBase.__init__(self,inp,pot)
        self.omega_a = inp.tani_omega_a

    def rho(self, qs):
        n_q = len(qs)
        m = self.mass
        om = self.omega_a
        beta = self.beta
        a = m*om/(2*hbar*sinh(beta*hbar*om))
        b = cosh(beta*hbar*om)

        rho0 = np.empty([n_q,n_q], dtype=complex)
        for i in range(n_q):
            rho0[i,:] = sqrt(a/pi)*np.exp(-a*b*(qs[i]**2 + qs**2)
                    + 2*a*qs[i]*qs)
        n_q = len(qs)
        dq = (qs[-1]-qs[0])/(n_q-1)
        return rho0*dq

#******************************************************************************
#******************************************************************************
                
class StateThermEq_aren(StateBase):
    """Thermal equilibrium isolated from the bath initial state"""
    def __init__(self,inp,pot):
        StateBase.__init__(self,inp,pot)
        self.pot0 = copy.deepcopy(pot)

    def rho(self,qs):
        n_q = len(qs)
        dq = (qs[-1]-qs[0])/(len(qs)-1)
        beta = self.beta
        print("Generating rhoq using potential that has a_ren: ", self.pot0.a_ren)
        ham, evals, evecs = dvr.ham_evals_evecs(qs, self.pot0, self.mass)
        n_EE = len(evals)
        # Calculating the partition function
        partfun = 0
        for n in range(n_EE):
            tmp = np.exp(-beta*(evals[n]-evals[0]))
            partfun += tmp
        # Generating ground state density
        # (in the bare Hamiltonian representation)
        rho = np.zeros([n_EE,n_EE],dtype=complex)
        count = 1
        for n in range(n_EE):
            tmp = np.exp(-beta*(evals[n]-evals[0]))/partfun
            rho[n,n] = tmp*np.conj(tmp)
            if n==0:
                norm = tmp
            else:
                if tmp>1e-5*norm:
                    count+=1
        print("The number of occupied states is ", count)

        # Converting to position representation
        rho = evecs@rho@evecs.T.conj()
        return rho

    def rho2(self,qs):
        # Equivalent, but not numerically identical to the function above
        n_q = len(qs)
        beta = self.beta
        ham, evals, evecs = dvr.ham_evals_evecs(qs, self.pot0, self.mass)
        n_EE = len(evals)
        # Calculating the partition function
        partfun = 0
        for n in range(n_EE):
            tmp = np.exp(-beta*(evals[n]-evals[0]))
            partfun += tmp
        rho = np.zeros([n_q,n_q],dtype=complex)
        count = 1
        for n in range(n_EE):
            tmp_coef = np.exp(-beta*(evals[n]-evals[0]))/partfun
            tmp = tmp_coef*np.outer(evecs[:,n],evecs[:,n])
            rho += tmp
            if n==0:
                norm = tmp_coef
            else:
                if tmp_coef>0.001*norm:
                    count+=1
        return rho

    def rho_q_kubo_q(self,qs):
        beta = self.beta
        n_q = len(qs)
        dq = (qs[-1]-qs[0])/(n_q-1)
        ham, evals, evecs = dvr.ham_evals_evecs(qs, self.pot0, self.mass)
        n_EE = len(evals)
        # Making q operator matrix in energy eigenvector representation
        q_mat = np.zeros([n_EE, n_EE])
        for i in range(n_EE):
            for j in range(n_EE):
                q_mat[i,j] = np.vdot(evecs[:,i], evecs[:,j]*qs)
        # Calculating the partition function
        partfun = 0
        for n in range(n_EE):
            tmp = np.exp(-beta*(evals[n]-evals[0]))
            partfun += tmp
        # Calculating the initial matrix in the EE representation
        rho = np.zeros([n_EE,n_EE],dtype=complex)
        count = 1
        for n in range(n_EE):
            for m in range(n_EE):
                if n==m:
                    rho[m,n] = q_mat[m,n]*np.exp(-beta*(evals[n]-evals[0]))/partfun
                else:
                    rho[m,n] = q_mat[m,n]*( (np.exp(-beta*(evals[n]-evals[0]))
                            -np.exp(-beta*(evals[m]-evals[0])))
                            /(beta*(evals[m]-evals[n])*partfun) )
                if n==0 and m==0:
                    norm = rho[m,n]
                else:
                    if tmp>1e-5*norm:
                        count+=1
        print("The number of occupied states in Kubo rho is ", count)
        # Converting to position representation
        #plt.contour(qs,qs,rho)
        #plt.show()
        #plt.clf()
        #fig = plt.figure()
        #ax = fig.gca(projection='3d')
        #X,Y = np.meshgrid(qs,qs)
        #surf = ax.plot_surface(X,Y, rho.real, antialiased=True, cmap=mpl.cm.coolwarm)
        #plt.show()
        #plt.clf()
        rho = evecs@rho@evecs.T.conj()

        #plt.contour(qs,qs,rho.real)
        #plt.show()
        #plt.clf()

        #plt.contour(qs,qs,rho.imag)
        #plt.show()
        #plt.clf()
        #plt.contour(qs,qs,np.abs(rho))
        #plt.show()
        #plt.clf()

        #fig = plt.figure()
        #ax = fig.gca(projection='3d')
        #X,Y = np.meshgrid(qs,qs)
        #surf = ax.plot_surface(X,Y, rho.real, antialiased=True, cmap=mpl.cm.coolwarm)
        #plt.show()
        #plt.clf()
        return rho

    def rho_q_kubo_q2(self,qs):
        beta = self.beta
        n_q = len(qs)
        dq = (qs[-1]-qs[0])/(n_q-1)
        ham, evals, evecs = dvr.ham_evals_evecs(qs, self.pot0, self.mass)
        n_EE = len(evals)
        # Making q operator matrix in energy eigenvector representation
        q_mat = np.zeros([n_EE, n_EE])
        for i in range(n_EE):
            for j in range(n_EE):
                q_mat[i,j] = np.vdot(evecs[:,i], evecs[:,j]*np.square(qs))
        # Calculating the partition function
        partfun = 0
        for n in range(n_EE):
            tmp = np.exp(-beta*(evals[n]-evals[0]))
            partfun += tmp
        # Calculating the initial matrix in the EE representation
        rho = np.zeros([n_EE,n_EE],dtype=complex)
        count = 1
        for n in range(n_EE):
            for m in range(n_EE):
                if n==m:
                    rho[m,n] = q_mat[m,n]*np.exp(-beta*(evals[n]-evals[0]))/partfun
                else:
                    rho[m,n] = q_mat[m,n]*( (np.exp(-beta*(evals[n]-evals[0]))
                            -np.exp(-beta*(evals[m]-evals[0])))
                            /(beta*(evals[m]-evals[n])*partfun) )
                if n==0 and m==0:
                    norm = rho[m,n]
                else:
                    if tmp>1e-5*norm:
                        count+=1
        print("The number of occupied states in Kubo rho is ", count)
        # Converting to position representation
        #plt.contour(qs,qs,rho)
        #plt.show()
        #plt.clf()
        #fig = plt.figure()
        #ax = fig.gca(projection='3d')
        #X,Y = np.meshgrid(qs,qs)
        #surf = ax.plot_surface(X,Y, rho.real, antialiased=True, cmap=mpl.cm.coolwarm)
        #plt.show()
        #plt.clf()
        rho = evecs@rho@evecs.T.conj()

        #plt.contour(qs,qs,rho.real)
        #plt.show()
        #plt.clf()

        #plt.contour(qs,qs,rho.imag)
        #plt.show()
        #plt.clf()
        #plt.contour(qs,qs,np.abs(rho))
        #plt.show()
        #plt.clf()

        #fig = plt.figure()
        #ax = fig.gca(projection='3d')
        #X,Y = np.meshgrid(qs,qs)
        #surf = ax.plot_surface(X,Y, rho.real, antialiased=True, cmap=mpl.cm.coolwarm)
        #plt.show()
        #plt.clf()
        return rho
#******************************************************************************
#******************************************************************************

class StateThermEq(StateThermEq_aren):
    """Thermal equilibrium of system isolated from the bath"""
    def __init__(self,inp,pot):
        StateThermEq_aren.__init__(self,inp,pot)
        self.pot0.a_ren = 0


#******************************************************************************
#******************************************************************************

class StateLoaded(StateBase):
    def __init__(self,inp,pot):
        StateBase.__init__(self,inp,pot)

    def rho(self, qs):
        filename = "rho0"
        limit = 1e-6
        if not os.path.isfile(filename):
            raise RuntimeError("Cannot find file " + filename
                                                    + "for initial state")
        with open(filename,'rb') as f:
            qs0, rho = pickle.load(f)
        # Checking compatibility of space grids
        if len(qs) != len(qs0):
            raise ValueError("Loaded rho0 has incompatible grid size")
        for i in range(len(qs)):
            if abs(qs[i]-qs0[i])>limit:
                raise ValueError("Loaded rho0 has incompatible grid")
        return rho

    def rho_EE(self, qs):
        raise NotImplementedError(
                "rho_EE not implemented for this initial state.")

    def wigner(self, qs):
        raise NotImplementedError(
                "wigner not implemented for this initial state.")

    def rho_q_kubo_q(self,qs):
        raise NotImplementedError(
                "rho_q_kubo_q not implemented for this initial state.")
    def rho_q_kubo_q2(self,qs):
        raise NotImplementedError(
                "rho_q_kubo_q2 not implemented for this initial state.")
#******************************************************************************
#******************************************************************************
# Analytic way of getting the wavepacket Wigner function
    #
    #n_q = len(qs)
    #n_p = len(qs)
    ## Creating ps
    #dp = dq*hbar/(2*self.var)
    #ps = np.zeros([len(qs)])
    #for i in range(n_p):
    #    #self.ps[i] = self.dp*(i-(self.n_p-1)/2)
    #    ps[i] = dp*(i+1-n_p/2) # Fortran compatibility
    ## Creating Wigner transform
    #varq = self.var
    #varp = hbar**2 / (4*self.var)
    #wig = np.empty([n_q,n_p])
    #for i in range(n_q):
    #    for j in range(n_p):
    #        wig[i,j] = (1/(pi*hbar)) * np.exp(
    #                -( (qs[i]-self.q0)**2 / (2*varq) )
    #                -( (ps[j]-self.p0)**2 / (2*varp) ))

# Analytic way of getting the Tanimura initial state Wigner function
    #
    #n_q = len(qs)
    #n_p = len(qs)
    ## Creating ps
    #dp = (dq*2*hbar*sqrt(cosh(self.beta*hbar*self.omega_a)**2 - 1)
    #        * self.mass*self.omega_a
    #            /(2*hbar*sinh(self.beta*hbar*self.omega_a)))
    #ps = np.zeros([len(qs)])
    #for i in range(n_p):
    #    #self.ps[i] = self.dp*(i-(self.n_p-1)/2)
    #    ps[i] = dp*(i+1-n_p/2) # Fortran compatibility
    ## Creating wig
    #m = self.mass
    #om = self.omega_a
    #beta = self.beta
    #a = m*om/(2*hbar*sinh(beta*hbar*om))
    #b = cosh(beta*hbar*om)
    #wig = np.empty([n_q,n_p])
    #for i in range(n_q):
    #    wig[i,:] = (1/(2*pi*hbar))*sqrt(2/(b+1))*np.exp(2*a*(1-b)*qs[i]**2
    #            - ps**2 / (2*a*(b+1)*hbar**2) )
