#!/usr/bin/env python
# File: simulations.py
"""Simulation objects.

Contains all simulation objects. Currently only
"""
from math import tan, sqrt
import sys
import numpy as np
import copy
import pickle

from heom.general import hbar, pi, formatflt, write_complex_pair, write_int_mat, write_real_mat, write_cplx_mat, get_flattening_coefficients
import heom.dvr as dvr
import heom.psd as psd
# Add the nbead jagged mode debye bath from my code
sys.path.append('/home/ach221/software/phd/heom/')
from debye_bath import Debye_bath 

#import pyheom
#******************************************************************************
#******************************************************************************

class SimBase(object):
    """ Base object for all simulation Sim objects"""
    def __init__(self, inp):
        # Reading in the input object
        self.sim_type = inp.sim_type
        self.bath_type = inp.bath_type
        self.sop_decomposition_type = inp.sop_decomposition_type
        self.init_state = inp.init_state
        self.pot_type = inp.pot_type
        self.truncation_type = inp.truncation_type
        # Cutoff coefficient per unit of time from the input file
        self.cutoff_coef = inp.cutoff_coef*inp.dt
        self.scaling_type = inp.scaling_type
        self.stability_threshold = inp.stability_threshold

        self.mass = inp.mass
        self.beta = inp.beta
        self.eta = inp.eta
        self.gamma = inp.gamma

        self.cL = inp.cL
        self.q0 = inp.q0
        self.q1 = inp.q1
        self.dq = inp.dq
        self.t0 = inp.t0
        # Taking into account that Kubo TCF can currently be calculated only
        # from a direct-product initial state
        if inp.heom_switch in ["kubo_qq","kubo_q2q2"]:
            self.t0 = 0
        self.t1 = inp.t1
        self.dt = inp.dt
        self.t_sample = inp.t_sample

        self.n_q = int((self.q1-self.q0)/self.dq)
        self.n_t = int((self.t1-self.t0)/self.dt)
        self.tau = int(self.t_sample/self.dt)
        if self.tau==0:
            self.tau = 1

        self.c0R = 0.5*hbar*self.eta*(self.gamma**2)/tan(
                                                0.5*self.beta*hbar*self.gamma)
        self.c0I = -0.5*hbar*self.eta*(self.gamma**2)
        self.c0 = self.c0R + 1j*self.c0I

        # Setting up the position grid
        self.ns = np.arange(0,self.cL+1)
        self.qs = np.zeros([self.n_q])
        for i in range(self.n_q):
            # + 1 is there just for consistency with Fortran
            self.qs[i] = self.dq*(i + 1 + int(self.q0/self.dq))


        # Setting up the time grid
        if self.t0<0 and self.t1>0:
            self.ts = np.zeros([self.n_t+1])
            n_neg = int(-self.t0/self.dt)
            for i in range(-n_neg, self.n_t+1-n_neg):
                self.ts[i+n_neg] = self.dt*i
            self.t_offset = n_neg%self.tau

        else:
            self.t_offset = 0
            self.ts = np.zeros([self.n_t+1])
            for i in range(self.n_t+1):
                self.ts[i] = self.dt*i + self.t0


    def write_position_grid(self, filename):
        with open(filename, "w+") as f:
            f.write(formatflt(0.0))
            for i in range(self.n_q):
                f.write(formatflt(self.qs[i]))
            f.write("\n")

    def save_checkpoint(self,filename,time):
        with open(filename,'wb') as f:
            pickle.dump((time, self.ados),f)

    def load_checkpoint(self,filename):
        with open(filename,'rb') as f:
            (time,self.ados) = pickle.load(f)
        self.t0 = time
        self.n_t = int((self.t1-self.t0)/self.dt)
        self.ts = np.zeros([self.n_t+1])
        for i in range(self.n_t+1):
            self.ts[i] = self.dt*i + self.t0

    def write_potential_grid(self, filename, pot):
        with open(filename, "w+") as f:
            for i in range(self.n_q):
                f.write(formatflt(self.qs[i]))
                f.write(formatflt(pot.calc(self.qs[i])))
                f.write("\n")

    def write_density(self, filename, time):
        raise SystemExit(
                "write_density not implemented for this simulation.")

    def write_q2(self, filename, time):
        raise SystemExit(
                "write_q2 not implemented for this simulation.")

    def check_stability(self, maxval):
        raise SystemExit(
                "Checking stability not implemented for this simulation.")

    def postmultiply_by_qs(self):
        raise SystemExit(
                "postmultiply_by_qs not implemented for this simulation.")

    def write_tcf_qq(self, filename_real, filename_imag, time):
        raise SystemExit(
                "write_tcf_qq not implemented for this simulation.")

    def write_tcf_q2q2(self, filename_real, filename_imag, time):
        raise SystemExit(
                "write_tcf_qq not implemented for this simulation.")

    def write_tcf_q4q4(self, filename_real, filename_imag, time):
        raise SystemExit(
                "write_tcf_qq not implemented for this simulation.")

    def write_ado_trend(self, filename, time):
        raise SystemExit(
                "write_ado_trend not implemented for this simulation.")

#******************************************************************************
#******************************************************************************
class SimQ(SimBase):
    """ Base object for all quantum HEOM simulation Sim objects"""
    def __init__(self,inp,pot,state):
        SimBase.__init__(self,inp)

        # Including -i/hbar in qs
        self.iqs = (-1j/hbar)*self.qs
        # Making Hamiltonian
        self.make_ham(pot)
        # Setting initial state
        self.rho0 = np.zeros([self.n_q, self.n_q], dtype=complex)
        self.rho0[:,:] = state.rho(self.qs)
        # Normalising the initial state
        print("Normalisation constant")
        print(abs(np.trace(self.rho0)))
        self.rho0 /= np.trace(self.rho0)

        # Q' matrix, i.e. for post-multiplying
        self.q_mat = np.zeros([self.n_q, self.n_q], dtype=complex)
        for j in range(self.n_q):
            self.q_mat[:,j] = self.qs[j]

    def make_ham(self, pot):
        self.ham = np.zeros([self.n_q, self.n_q])
        for i in range(self.n_q):
            for j in range(self.n_q):
                if i==j:
                   self.ham[i,j] = ( pot.calc(self.qs[i]) 
                           + hbar**2 *  pi**2 / (6*self.mass*self.dq**2) )
                else:
                   self.ham[i,j] = ( hbar**2 * (-1)**(i-j)
                       / (self.mass * self.dq**2 * (i-j)**2) )
        # Below should replace the lines above
        #self.ham = dvr.hamiltonian(self.qs,pot,self.mass)
        self.iham = (-1j/hbar)*self.ham

    def update_pot(self, pot):
        self.make_ham(pot)

#******************************************************************************
#******************************************************************************

class SimQHEOM(SimQ):
    """ Quantum HEOM simulation without Matsubara terms"""
    def __init__(self, inp, pot, state):
        SimQ.__init__(self,inp, pot, state)

        self.ados = np.zeros([self.cL+1, self.n_q, self.n_q], dtype=complex)
        self.ados[0,:,:] = np.copy(self.rho0)

        self.k1 = np.zeros_like(self.ados)
        self.k2 = np.zeros_like(self.ados)
        self.k3 = np.zeros_like(self.ados)
        self.k4 = np.zeros_like(self.ados)

        self.q_rho = np.zeros_like(self.ados)
        self.rho_q = np.zeros_like(self.ados)
        # self.comm = np.zeros_like(self.ados)
        # self.anticomm = np.zeros_like(self.ados)
        # Scaling variables
        self.coef_up = np.ones([self.cL+1])
        self.coef_down = np.ones([self.cL+1])
        if self.scaling_type=="shi":
            for nk in range(self.cL+1):
                self.coef_up[nk] = sqrt((nk+1)*abs(self.c0))
                if nk!=0:
                    self.coef_down[nk] =  sqrt(nk/abs(self.c0))/nk

    def get_kubo_q_state(self, state):
        self.ados = np.zeros([self.cL+1, self.n_q, self.n_q], dtype=complex)
        self.ados[0,:,:] = state.rho_q_kubo_q(self.qs)
        self.k1 = np.zeros_like(self.ados)
        self.k2 = np.zeros_like(self.ados)
        self.k3 = np.zeros_like(self.ados)
        self.k4 = np.zeros_like(self.ados)
    def get_kubo_q2_state(self, state):
        self.ados = np.zeros([self.cL+1, self.n_q, self.n_q], dtype=complex)
        self.ados[0,:,:] = state.rho_q_kubo_q2(self.qs)
        self.k1 = np.zeros_like(self.ados)
        self.k2 = np.zeros_like(self.ados)
        self.k3 = np.zeros_like(self.ados)
        self.k4 = np.zeros_like(self.ados)

    def propagate(self, pot, dt):
        # Disguised Runge-Kutta
        # k1 = dt*derivative(rho)
        # k2 = dt*derivative(rho+k1/2)
        # k3 = dt*derivative(rho+k2/2)
        # k4 = dt*derivative(rho+k3)
        # rho = rho + (k1+k4)/6 + (k2+k3)/3

        self.derivative(self.ados, self.k1)
        self.k1 = (dt/2)*self.k1
        self.derivative(self.ados+self.k1, self.k2)
        self.k2 = (dt/2)*self.k2
        self.derivative(self.ados+self.k2, self.k3)
        self.k3 = dt*self.k3
        self.derivative(self.ados+self.k3, self.k4)
        self.ados += (dt/6)*self.k4 + (2/3)*self.k2 + (self.k3+self.k1)/3


    def derivative(self, rho, drdt):
        # Using rho for ados for compatibility reason
        if self.bath_type=="none":
            n=0 # Zeroth tier ADO
            drdt[n,:,:] = self.iham@rho[n,:,:] - rho[n,:,:]@self.iham

        elif self.bath_type=="debye":
            ## Naive implementation
            ##for n in range(self.cL+1):
            ##    for i in range(self.n_q):
            ##        for j in range(self.n_q):
            ##            self.q_rho[n,i,j] = rho[n,i,j]*self.qs[i]
            ##            self.rho_q[n,i,j] = rho[n,i,j]*self.qs[j]

            ## Elegant, but not as fast implementation (even with optimize)
            ##self.q_rho = np.einsum('...ij,i->...ij',rho,self.qs)
            ##self.rho_q = np.einsum('...ij,j->...ij',rho,self.qs)

            # NB: q_rho and rho_q alreadu include -i/hbar
            self.q_rho = np.swapaxes(
                    np.multiply(np.swapaxes(rho,1,2), self.iqs),1,2)
            self.rho_q = np.multiply(rho, self.iqs)
            #*****************************************************************
            # VERSION 1: Less efficient, but more explicit
            #self.comm[:,:,:] = (-1j/hbar)*(self.q_rho - self.rho_q)
            #self.anticomm[:,:,:] = self.q_rho + self.rho_q
            ## NB: iham includes -i/hbar
            #n=0 # Zeroth tier ADO
            #drdt[n,:,:] = (self.iham@rho[n,:,:] - rho[n,:,:]@self.iham
            #        + self.comm[n+1,:,:])
            #n=self.cL # Highest tier ADO
            #drdt[n,:,:] = (
            #          self.iham@rho[n,:,:] - rho[n,:,:]@self.iham
            #        - n*self.gamma*rho[n,:,:]
            #        + (n*self.c0R)*self.comm[n-1,:,:]
            #        - (n*self.eta*self.gamma/2)*self.anticomm[n-1,:,:])
            #for n in range(1, self.cL):
            #    drdt[n,:,:] = (
            #            self.iham@rho[n,:,:] - rho[n,:,:]@self.iham
            #            - n*self.gamma*rho[n,:,:]
            #            + self.comm[n+1,:,:]
            #            + (n*self.c0R)*self.comm[n-1,:,:]
            #            - (n*self.eta*self.gamma/2)*self.anticomm[n-1,:,:])
            #*****************************************************************
            # VERSION 2:
            # NB: iham includes -i/hbar
            if self.cL==0:
                n=0 # Zeroth tier ADO
                drdt[n,:,:] = self.iham@rho[n,:,:] - rho[n,:,:]@self.iham
            elif self.cL>0:
                n=0 # Zeroth tier ADO
                drdt[n,:,:] = (self.iham@rho[n,:,:] - rho[n,:,:]@self.iham
                            + self.coef_up[n]*(self.q_rho[n+1,:,:]
                            - self.rho_q[n+1,:,:]))
                n=self.cL # Highest tier ADO
                drdt[n,:,:] = (
                          self.iham@rho[n,:,:] - rho[n,:,:]@self.iham
                        - n*self.gamma*rho[n,:,:]
                        + (self.coef_down[n]*n*self.c0)*self.q_rho[n-1,:,:]
                        - (self.coef_down[n]*n*self.c0.conjugate())*self.rho_q[n-1,:,:]
                        )
            if self.cL>1:
                for n in range(1, self.cL):
                    drdt[n,:,:] = (
                            self.iham@rho[n,:,:] - rho[n,:,:]@self.iham
                            - n*self.gamma*rho[n,:,:]
                            + self.coef_up[n]
                                *(self.q_rho[n+1,:,:] - self.rho_q[n+1,:,:])
                            + (self.coef_down[n]*n*self.c0)*self.q_rho[n-1,:,:]
                            - (self.coef_down[n]*n*self.c0.conjugate())
                                                            *self.rho_q[n-1,:,:]
                            )
        else:
            raise SystemExit("Chosen bath type not supported for this sim")

    def check_stability(self, maxval):
        outcome = True
        for i in range(self.n_q):
            if abs(self.ados[0,i,i])>maxval:
                outcome = False
                break
        return outcome

    def postmultiply_by_qs(self):
        #self.ados[0,:,:] *= self.q_mat
        # Multiplying all ados
        for n in range(self.cL):
            self.ados[n,:,:] *= self.q_mat


    def write_tcf_qq(self, filename_real, filename_imag, time):
        write_complex_pair(filename_real, filename_imag, time,
                        np.trace(self.ados[0,:,:]*self.q_mat))

    def write_tcf_q2q2(self, filename_real, filename_imag, time):
        write_complex_pair(filename_real, filename_imag, time,
                        np.trace(self.ados[0,:,:]*self.q_mat*self.q_mat))

    def write_tcf_q4q4(self, filename_real, filename_imag, time):
        write_complex_pair(filename_real, filename_imag, time,
        np.trace(self.ados[0,:,:]*self.q_mat*self.q_mat*self.q_mat*self.q_mat)
            )

    def write_density(self, filename, time):
        with open(filename, "a") as f:
            f.write(formatflt(time))
            for i in range(self.n_q):
                f.write(formatflt(np.real(self.ados[0,i,i])/self.dq))
            f.write("\n")

    def write_ado_trend(self, filename, time):
        with open(filename, "a") as f:
            f.write(formatflt(time))
            for i in range(self.cL+1):
                f.write(formatflt(np.amax(abs(self.ados[i,:,:]))))
            f.write("\n")

#******************************************************************************
#******************************************************************************

class SimC(SimBase):
    """ Base object for all classical HEOM simulation Sim objects"""
    def __init__(self,inp,pot,state):
        SimBase.__init__(self,inp)
        self.n_p = self.n_q
        print("getting initial state")
        self.wig0, self.ps, self.dp = state.wigner(self.qs,self.dq)
        #self.dp = state.dp(self.dq)
        #self.ps = np.zeros([self.n_p])
        #for i in range(self.n_p):
        #    #self.ps[i] = self.dp*(i-(self.n_p-1)/2)
        #    self.ps[i] = self.dp*(i+1-self.n_p/2) # Fortran compatibility
        #self.wig0 = state.wigner(self.qs, self.ps)

        # Normalising the initial state
        print("Normalisation constant")
        print(self.dq*self.dp*np.sum(self.wig0))
        self.wig0 /= self.dq*self.dp*np.sum(self.wig0)
        self.dq_mat = dvr.derivative1(self.qs)
        self.dp_mat = dvr.derivative1(self.ps)

#******************************************************************************
#******************************************************************************

class SimCHEOM(SimC):
    """ Classical HEOM simulation without Matsubara terms"""
    def __init__(self, inp, pot, state):
        SimC.__init__(self,inp, pot, state)
        self.ados = np.zeros([self.cL+1, self.n_q, self.n_p])
        self.ados[0,:,:] = np.copy(self.wig0)

        self.k1 = np.zeros_like(self.ados)
        self.k2 = np.zeros_like(self.ados)
        self.k3 = np.zeros_like(self.ados)
        self.k4 = np.zeros_like(self.ados)

        self.ddq = np.zeros_like(self.ados)
        self.ddp = np.zeros_like(self.ados)

        self.p_mat = np.tile(self.ps, (self.cL+1, self.n_q, 1))
        self.q_mat = np.zeros_like(self.ados)
        for i in range(self.n_q):
            self.q_mat[:,i,:] = self.qs[i]
        self.make_pot_mat(pot)
        self.coef_up = np.ones([self.cL+1])
        self.coef_down = np.ones([self.cL+1])
        if self.scaling_type=="shi":
            for nk in range(self.cL+1):
                self.coef_up[nk] = sqrt((nk+1)*abs(self.c0))
                if nk!=0:
                    self.coef_down[nk] =  sqrt(nk/abs(self.c0))/nk

    def make_pot_mat(self,pot):
        self.pot_mat = np.zeros([self.cL+1, self.n_q, self.n_p])
        for n in range(self.cL+1):
            for i in range(self.n_q):
                self.pot_mat[n,i,:] = pot.diff(self.qs[i])

    def update_pot(self,pot):
        self.make_pot_mat(pot)

    def propagate(self, pot, dt):
        # Disguised Runge-Kutta
        # k1 = dt*derivative(rho)
        # k2 = dt*derivative(rho+k1/2)
        # k3 = dt*derivative(rho+k2/2)
        # k4 = dt*derivative(rho+k3)
        # rho = rho + (k1+k4)/6 + (k2+k3)/3
        self.k1[:,:,:] = 0.0 
        self.k2[:,:,:] = 0.0 
        self.k3[:,:,:] = 0.0 
        self.k4[:,:,:] = 0.0 

        self.derivative(self.ados, self.k1)
        self.k1 = (dt/2)*self.k1
        self.derivative(self.ados+self.k1, self.k2)
        self.k2 = (dt/2)*self.k2
        self.derivative(self.ados+self.k2, self.k3)
        self.k3 = dt*self.k3
        self.derivative(self.ados+self.k3, self.k4)
        self.ados += (dt/6)*self.k4 + (2/3)*self.k2 + (self.k3+self.k1)/3

    def derivative(self, wig, drdt):
        # Calculating derivative matrices
        for n in range(self.cL+1):
            for i in range(self.n_q):
                self.ddp[n,i,:] = self.dp_mat@wig[n,i,:]
            for j in range(self.n_p):
                self.ddq[n,:,j] = self.dq_mat@wig[n,:,j]
        # This is MINUS the Liouvillian
        lvln = (-1/self.mass)*self.p_mat*self.ddq + self.ddp*self.pot_mat

        if self.bath_type=="none":
            n=0 # Zeroth tier ADO
            drdt[0,:,:] = lvln[0,:,:]

        elif self.bath_type in {"debye","debye_correction3_cl"}:
            if self.bath_type == "debye":
                theta_f = ((self.eta*self.gamma/self.beta)
                        *(self.ddp-self.beta*self.gamma*self.q_mat*wig))
            elif self.bath_type == "debye_correction3_cl":
                theta_f = ((self.eta*self.gamma**2/self.mass)*self.p_mat*wig
                    + (self.eta*self.gamma**2/self.beta)*self.ddp)
            if self.cL ==0:
                n=0 # Zeroth tier ADO
                drdt[n,:,:] = lvln[n,:,:]
            elif self.cL>0:
                # Zeroth tier ADO
                n=0
                drdt[n,:,:] = lvln[n,:,:] + self.coef_up[n]*self.ddp[n+1,:,:]
                # Highest tier ADO
                n=self.cL
                drdt[n,:,:] = (
                        + lvln[n,:,:]
                        - n*self.gamma*wig[n,:,:]
                        + self.coef_down[n]*n*theta_f[n-1,:,:]
                        )
            if self.cL>1:
                for n in range(1, self.cL):
                    drdt[n,:,:] = (
                            lvln[n,:,:]
                            - n*self.gamma*wig[n,:,:]
                            + self.coef_up[n]*self.ddp[n+1,:,:]
                            + self.coef_down[n]*n*theta_f[n-1,:,:]
                            )
        else:
            raise SystemExit("Chosen bath type not supported for this sim")

    def write_density(self, filename, time):
        with open(filename, "a") as f:
            f.write(formatflt(time))
            data = self.dp*np.sum(self.ados[0,:,:],axis=1)
            for i in range(self.n_q):
                f.write(formatflt(data[i]))
            f.write("\n")

    def check_stability(self, maxval):
        outcome = True
        if np.amax(self.dp*np.sum(self.ados[0,:,:],axis=1))>maxval:
                outcome = False
        return outcome

    def write_ado_trend(self, filename, time):
        with open(filename, "a") as f:
            f.write(formatflt(time))
            for i in range(self.cL+1):
                f.write(formatflt(np.amax(abs(self.ados[i,:,:]))))
            f.write("\n")

    def postmultiply_by_qs(self):
        #self.ados[0,:,:] *= self.q_mat[0,:,:]
        for n in range(self.cL):
            # This may need to be checked
            self.ados[n,:,:] *= self.q_mat[0,:,:]

    def write_tcf_qq(self, filename_real, filename_imag, time):
        write_complex_pair(filename_real, filename_imag, time,
            self.dq*self.dp*np.sum(self.ados[0,:,:]*self.q_mat[0,:,:]))

    def write_tcf_q2q2(self, filename_real, filename_imag, time):
        write_complex_pair(filename_real, filename_imag, time,
            self.dq*self.dp*np.sum(self.ados[0,:,:]*self.q_mat[0,:,:]*self.q_mat[0,:,:])
            )

    def write_tcf_q4q4(self, filename_real, filename_imag, time):
        write_complex_pair(filename_real, filename_imag, time,
            self.dq*self.dp*np.sum(self.ados[0,:,:]*self.q_mat[0,:,:]
                    *self.q_mat[0,:,:]*self.q_mat[0,:,:]*self.q_mat[0,:,:])
            )

#******************************************************************************
#******************************************************************************

class Tensor(dict):
    """Dictionary-based object for storing ADOs"""
    # Dictionary based tensor object
    # no need for overriding __init__
    def insert(self, key, rho):
        # Inserts ADO, either adding or creating a new item
        tup = tuple(key)
        if tup in self:
            self[tup] += rho
        else:
            self[tup] = np.copy(rho)

    def remove(self, key):
        # Removes item if present and does not complain if not
        self.pop(tup, None) 
        tup = tuple(key)

    def add(self, *args):
        for other in args:
            for n, rho in other.items():
                if n in self:
                    self[n] += rho
                else:
                    self[n] = np.copy(rho)
        return self

    def __add__(self, other):
        dic = copy.deepcopy(self)
        for n, rho in other.items():
            if n in dic:
                dic[n] += rho
            else:
                dic[n] = np.copy(rho)
        return dic

    def times(self, coef):
        for n in self:
            self[n] *= coef
        return self

    def prune(self, cutoff):
        n_list = []
        for n, rho in self.items():
            if np.amax(np.abs(rho))<cutoff:
                n_list.append(n)
        for n in n_list:
            if sum(n)!=0: # Never delete physical density matrix
                del self[n]

#******************************************************************************
#******************************************************************************

class SimMats(object):
    """Matsubara terms object

       This object carries variables and methods for HEOM codes that
       include the Matsubara terms.
       """
    def __init__(self,inp):
        # Assumes being called with self containing necessary variables
        self.pruning_type = inp.pruning_type
        self.cK = inp.cK
        self.Ks = np.arange(start=0,stop=self.cK+1,step=1, dtype=int)
        self.cL = self.cL

        self.gammas = self.make_gammas()
        self.make_cs() # makes self.cs, self.low_T_coef_R/I
        self.n0 = tuple(np.zeros(self.cK+1, dtype=int))
        # Scaling variables
        self.coef_up = np.ones([self.cL+1,self.cK+1])
        self.coef_down = np.ones([self.cL+1,self.cK+1])
        if self.scaling_type=="shi":
            for nk in range(self.cL+1):
                for k in range(self.cK+1):
                    self.coef_up[nk,k] = sqrt((nk+1)*abs(self.cs[k]))
                    if nk!=0:
                        self.coef_down[nk,k] =  sqrt(nk/abs(self.cs[k]))/nk

        self.qs2 = self.qs**2

    def make_gammas(self):
        gammas = []
        if "debye" in self.bath_type:
            gammas.append(self.gamma) #the zeroth gamma is the system one
            if self.sop_decomposition_type=="matsubara":
                for k in range(1,self.cK+1):
                    gammas.append(2*pi*k/(self.beta*hbar))
            elif "pade" in self.sop_decomposition_type:
                psd.scheme = self.sop_decomposition_type
                self.psd_xi, self.psd_eta = psd.get_xi_eta(self.cK)
                for k in range(1,self.cK+1):
                    gammas.append(self.psd_xi[k-1]/(self.beta*hbar))
                    # xi_k = beta*hbar*gamma_k
            elif self.sop_decomposition_type=="nbead": #NOTE- nbead BCF implementaion, sines go here
                bathmode = ['nbead','matsubara'][0]
                bathobj = Debye_bath(self.eta*self.gamma,self.gamma,self.beta,hbar,self.cK,bathmode)
                _,gammas = bathobj.get_coeffs() 
                gammas = gammas.tolist()
                # sys.exit()
        elif "white" in self.bath_type:
            if "pade" in self.sop_decomposition_type:
                psd.scheme = self.sop_decomposition_type
                self.psd_xi, self.psd_eta = psd.get_xi_eta(self.cK+1)
                for k in range(len(self.psd_xi)):
                    gammas.append(self.psd_xi[k]/(self.beta*hbar))
                    # xi_k = beta*hbar*gamma_k
            else:
                raise SystemExit("Unknown decomposition for white bath")
        elif "none"==self.bath_type:
            pass

        else:
            raise SystemExit("This bath has not yet been implemented")
        return tuple(gammas)

    def make_cs(self):
        self.cs = []
        #********************************************************
        # Debye
        #********************************************************
        if "debye" in self.bath_type:
            #********************************************************
            # Matsubara
            #********************************************************
            if self.sop_decomposition_type=="matsubara":
                self.cs.append(self.c0)
                for k in range(1, self.cK+1):
                    self.cs.append(
                            (2*self.eta*(self.gamma**2)/self.beta)
                            * self.gammas[k]/(self.gammas[k]**2-self.gamma**2)
                            )
                # Low temperature correction variables
                self.low_T_coef_R = self.eta/(self.beta*hbar**2)
                self.low_T_coef_R += (-1/(hbar**2))*self.c0R/self.gamma
                for i in range(1,self.cK+1):
                    self.low_T_coef_R += (-1/(hbar**2))*self.cs[i]/self.gammas[i]
                self.low_T_coef_I = (-1j/(hbar**2))*self.c0I/self.gamma
            #********************************************************
            # Nbead jagged mode
            #********************************************************
            elif self.sop_decomposition_type=="nbead":  #NOTE- nbead BCF implementaion, Cks go here
                bathmode = ['nbead','matsubara'][0]
                bathobj = Debye_bath(self.eta*self.gamma,self.gamma,self.beta,hbar,self.cK,bathmode)
                self.cs,_ = bathobj.get_coeffs() 
                self.cs = self.cs.tolist()
                self.low_T_coef_R = 0
                self.low_T_coef_I = 0
                #still not sure what the lowT correction will be (is this the xi thing?)
            #********************************************************
            # Pade
            #********************************************************
            elif "pade" in self.sop_decomposition_type:
                self.c0R = 1
                for k in range(1, self.cK+1):
                    self.c0R -= self.psd_eta[k-1]*2*(self.gamma**2) /(self.gammas[k]**2-self.gamma**2)
                self.c0R *= self.eta*self.gamma/self.beta
                self.c0 = self.c0R + 1j*self.c0I
                self.cs.append(self.c0)
                for k in range(1, self.cK+1):
                    self.cs.append(
                            (self.psd_eta[k-1]*2*self.eta*(self.gamma**2)/self.beta)
                            * self.gammas[k]/(self.gammas[k]**2-self.gamma**2)
                            )
                # Low temperature correction variables
                self.low_T_coef_R = 0
                # ^ If the Matsubara definition was used, it would give zero
                # anyway (within numerical precision)
                self.low_T_coef_I = 0
        #********************************************************
        # White
        #********************************************************
        elif "white" in self.bath_type:
            if "pade" in self.sop_decomposition_type:
                for k in range(self.cK+1):
                    self.cs.append(
                            -self.psd_eta[k]*2*self.eta*self.gammas[k]/self.beta
                            )
                # Low temperature correction variables
                # In this case really the delta-function component of alpha,
                # the bath collective variable TCF.

                # This is c_delta from Ikeda Tanimura 2019 JCTC
                self.low_T_coef_R = 1
                for k in range(self.cK+1):
                    self.low_T_coef_R += 2*self.psd_eta[k]
                self.low_T_coef_R *= self.eta/(self.beta)
                # This is r_delta
                self.low_T_coef_I = self.eta/(2*self.mass)
            else:
                print("Unknown decomposition for white bath")
                exit()
        elif "none" == self.bath_type:
            # Just a placeholder
            self.cs = np.zeros(self.cK+1)
        else:
            print("This bath has not yet been implemented")
        self.cs = tuple(self.cs)

    def propagate(self, pot, dt):
        # Disguised Runge-Kutta
        # k1 = dt*derivative(rho)
        # k2 = dt*derivative(rho+k1/2)
        # k3 = dt*derivative(rho+k2/2)
        # k4 = dt*derivative(rho+k3)
        # rho = rho + (k1+k4)/6 + (k2+k3)/3

        if self.pruning_type == "none":
            self.derivative(self.ados, dt/2, drdt=self.k1)
            self.derivative(self.ados+self.k1, dt/2, drdt=self.k2)
            self.derivative(self.ados+self.k2, dt, drdt=self.k3)
            self.derivative(self.ados+self.k3, dt/6, drdt=self.k4)
            self.k1.times(1/3)
            self.k2.times(2/3)
            self.k3.times(1/3)
            self.ados.add(self.k1,self.k2,self.k3,self.k4)
        elif self.pruning_type == "pruning":
            self.derivative(self.ados, dt/2, drdt=self.k1)
            self.derivative(self.ados+self.k1, dt/2, drdt=self.k2)
            self.derivative(self.ados+self.k2, dt, drdt=self.k3)
            self.derivative(self.ados+self.k3, dt/6, drdt=self.k4)
            self.k1.times(1/3)
            self.k2.times(2/3)
            self.k3.times(1/3)
            self.ados.add(self.k1,self.k2,self.k3,self.k4)
            self.ados.prune(self.cutoff)
        elif self.pruning_type == "radical":
            self.derivative(self.ados, dt/2, drdt=self.k1)
            self.k1.prune(self.cutoff)
            self.derivative(self.ados+self.k1, dt/2, drdt=self.k2)
            self.k2.prune(self.cutoff)
            self.derivative(self.ados+self.k2, dt, drdt=self.k3)
            self.k3.prune(self.cutoff)
            self.derivative(self.ados+self.k3, dt/6, drdt=self.k4)
            self.k4.prune(self.cutoff)
            self.k1.times(1/3)
            self.k2.times(2/3)
            self.k3.times(1/3)
            self.ados.add(self.k1,self.k2,self.k3,self.k4)
            self.ados.prune(self.cutoff)
        else:
            raise SystemExit("Unknown pruning method selected")


class SimQHEOM_Mats(SimQ,SimMats):
    """ Quantum HEOM simulation with Matsubara terms"""
    def __init__(self, inp, pot, state):
        SimQ.__init__(self,inp, pot, state)
        SimMats.__init__(self,inp)
        # Setting up the initial state
        self.cutoff = np.amax(np.abs(self.rho0))*self.cutoff_coef
        self.ados = Tensor()
        self.ados.insert(self.n0,self.rho0)

        # Work variables
        self.k1 = Tensor()
        self.k2 = Tensor()
        self.k3 = Tensor()
        self.k4 = Tensor()
        self.rho_q = np.zeros([self.n_q, self.n_q], dtype=complex)
        self.q_rho = np.zeros([self.n_q, self.n_q], dtype=complex)
        self.dq_mat = dvr.derivative1(self.qs)

    def get_kubo_q_state(self, state):
        new_rho = state.rho_q_kubo_q(self.qs)
        self.ados = Tensor()
        self.ados.insert(self.n0,new_rho)
        self.k1 = Tensor()
        self.k2 = Tensor()
        self.k3 = Tensor()
        self.k4 = Tensor()
    def get_kubo_q2_state(self, state):
        new_rho = state.rho_q_kubo_q2(self.qs)
        self.ados = Tensor()
        self.ados.insert(self.n0,new_rho)
        self.k1 = Tensor()
        self.k2 = Tensor()
        self.k3 = Tensor()
        self.k4 = Tensor()


    def derivative(self, ados, coef, drdt):
        drdt.clear()
        for n, rho_raw in ados.items():
            rho = coef*rho_raw
            if self.bath_type =="none":
                # Terms: Liouvillian + friction i.e. n = n
                drdt.insert(n,
                            self.iham@rho - rho@self.iham
                        )
                # NB: iham includes -i/hbar
            elif self.bath_type in {"debye","debye_correction","debye_correction2","white"}:
                if self.bath_type=="debye":
                    # Terms: Liouvillian + friction i.e. n = n
                    drdt.insert(n,
                                self.iham@rho - rho@self.iham
                                - np.dot(n,self.gammas)*rho
                            )
                    # NB: iham includes -i/hbar
                elif self.bath_type in {"debye_correction"}:
                    # Terms: Liouvillian + friction + low T correction
                    drdt.insert(n,
                                self.iham@rho - rho@self.iham
                                - np.dot(n,self.gammas)*rho
                                - self.low_T_coef_R*(
                                    np.swapaxes( # qq_rho -2*q_rho_q + rho_qq
                                        np.multiply(np.swapaxes(rho,0,1), self.qs2),0,1) 
                                    - 2*np.multiply(np.swapaxes(np.multiply(np.swapaxes(rho,0,1), self.qs),0,1), self.qs) 
                                    + np.multiply(rho, self.qs2) 
                                    )
                                )
                    # NB: iham includes -i/hbar
                elif self.bath_type=="debye_correction2":
                    # Terms: Liouvillian + friction + alternative low T correction
                    drdt.insert(n,
                                self.iham@rho - rho@self.iham
                                - np.dot(n,self.gammas)*rho
                                # -(R+I)*qq_rho -(R-I)*rho_qq -R*(-2*q_rho_q)
                                - (self.low_T_coef_R+self.low_T_coef_I)*np.swapaxes( np.multiply(np.swapaxes(rho,0,1), self.qs2),0,1)
                                - (self.low_T_coef_R-self.low_T_coef_I)*np.multiply(rho, self.qs2)
                                + 2*self.low_T_coef_R*np.multiply(np.swapaxes( np.multiply(np.swapaxes(rho,0,1), self.qs),0,1), self.qs)
                            )
                    # NB: iham includes -i/hbar
                elif self.bath_type in {"white"}:
                    # Note that the low_T_coef variables are used for coefficients
                    # of delta function terms in the bath TCF decomposition
                    drdt.insert(n,
                                self.iham@rho - rho@self.iham
                                - np.dot(n,self.gammas)*rho
                                - self.low_T_coef_R*(
                                    np.swapaxes( # qq_rho -2*q_rho_q + rho_qq
                                        np.multiply(np.swapaxes(rho,0,1), self.qs2),0,1) 
                                    - 2*np.multiply(np.swapaxes(np.multiply(np.swapaxes(rho,0,1), self.qs),0,1), self.qs) 
                                    + np.multiply(rho, self.qs2) 
                                    )
                                - self.low_T_coef_I*(
                                                self.q_mat@self.dq_mat@rho
                                                +self.q_mat@rho@self.dq_mat
                                                -self.dq_mat@rho@self.q_mat
                                                -rho@self.dq_mat@self.q_mat
                                                #self.q_mat@self.dq_mat@rho
                                                #-self.q_mat@rho@self.dq_mat
                                                #-self.dq_mat@rho@self.q_mat
                                                #+rho@self.dq_mat@self.q_mat
                                    )
                            )
                else:
                    raise SystemExit("Chosen bath type not supported for this sim")

                # Calculating the "Commutator and anticommutator"
                # NB: iqs = -i/hbar * qs
                self.q_rho = np.swapaxes(
                        np.multiply(np.swapaxes(rho,0,1), self.iqs),0,1)
                self.rho_q = np.multiply(rho, self.iqs)
                comm = self.q_rho-self.rho_q
                # rho_+k
                if sum(n)>0:
                    for k in range(self.cK+1):
                        if n[k]!=0:
                            n_temp = list(n)
                            n_temp[k] -= 1
                            drdt.insert(n_temp,self.coef_up[n_temp[k],k]*comm)
                # rho_-k
                if (sum(n)<self.cL) and (self.eta != 0):
                    # k > 0
                    for k in range(1,self.cK+1):
                        if n[k] != self.cL:
                            n_temp = list(n)
                            n_temp[k] += 1
                            drdt.insert(n_temp,
                                    self.coef_down[n_temp[k],k]
                                    * n_temp[k]*self.cs[k]*comm
                                    )
                    # k = 0
                    if n[0] != self.cL:
                        n_temp = list(n)
                        n_temp[0] += 1
                        drdt.insert(n_temp,
                                (self.coef_down[n_temp[0],0] * n_temp[0])
                                * (self.cs[0]*self.q_rho
                                    - self.cs[0].conjugate()*self.rho_q)
                                )
            else:
                raise SystemExit("Chosen bath type not supported for this sim")


    def write_density(self, filename, time):
        with open(filename, "a") as f:
            f.write(formatflt(time))
            for i in range(self.n_q):
                f.write(formatflt(np.real(self.ados[self.n0][i,i])/self.dq))
            f.write("\n")

    def postmultiply_by_qs(self):
        #self.ados[self.n0] *= self.q_mat
        # Multiplying all ADOs by the q matrix
        for key in self.ados:
            self.ados[key] *= self.q_mat

    def write_tcf_qq(self, filename_real, filename_imag, time):
        write_complex_pair(filename_real, filename_imag, time,
            np.trace(self.ados[self.n0]*self.q_mat))

    def write_tcf_q2q2(self, filename_real, filename_imag, time):
        write_complex_pair(filename_real, filename_imag, time,
            np.trace(self.ados[self.n0]*self.q_mat*self.q_mat))

    def write_tcf_q4q4(self, filename_real, filename_imag, time):
        write_complex_pair(filename_real, filename_imag, time,
            np.trace(self.ados[self.n0]*self.q_mat*self.q_mat*self.q_mat*self.q_mat)
            )

    def check_stability(self, maxval):
        outcome = True
        for i in range(self.n_q):
            if abs(self.ados[self.n0][i,i])>maxval:
                outcome = False
                break
        return outcome

    def write_ado_trend(self, filename, time):
        outcome = np.zeros([self.cL+1])
        for n, rho in self.ados.items():
            outcome[sum(n)] += np.amax(abs(self.ados[self.n0]))
        with open(filename, "a") as f:
            f.write(formatflt(time))
            for i in range(self.cL+1):
                f.write(formatflt(outcome[i]))
            f.write("\n")


#******************************************************************************
#******************************************************************************


class SimCHEOM_Mats(SimC,SimMats):
    """ Classical HEOM simulation with Matsubara terms"""
    def __init__(self, inp, pot, state):
        SimC.__init__(self,inp, pot, state)
        SimMats.__init__(self,inp)

        # Setting up the initial state
        self.ados = Tensor()
        self.ados.insert(self.n0,self.wig0)
        self.cutoff = np.amax(np.abs(self.wig0))*self.cutoff_coef

        # Work variables
        self.k1 = Tensor()
        self.k2 = Tensor()
        self.k3 = Tensor()
        self.k4 = Tensor()
        self.ddq = np.zeros([self.n_q, self.n_p])
        self.ddp = np.zeros([self.n_q, self.n_p])

        self.p_mat = np.tile(self.ps, (self.n_q, 1))
        self.q_mat = np.zeros([self.n_q, self.n_p])
        for i in range(self.n_q):
            self.q_mat[i,:] = self.qs[i]

        self.make_pot_mat(pot)

    def make_pot_mat(self,pot):
        self.pot_mat = np.zeros([self.n_q, self.n_p])
        for i in range(self.n_q):
            self.pot_mat[i,:] = pot.diff(self.qs[i])

    def update_pot(self,pot):
        self.make_pot_mat(pot)


    def derivative(self, ados, coef, drdt):
        drdt.clear()
        for n, wig_raw in ados.items():
            wig = coef*wig_raw
            # Calculating derivative matrices
            for i in range(self.n_q):
                self.ddp[i,:] = self.dp_mat@wig[i,:]
            for j in range(self.n_p):
                self.ddq[:,j] = self.dq_mat@wig[:,j]
            # This is MINUS the Liouvillian
            lvln = (-1/self.mass)*self.p_mat*self.ddq + self.ddp*self.pot_mat
            # Tanimura theta
            #theta_f = ((self.eta*self.gamma/self.mass)*self.p_mat*wig
            #            + (self.eta*self.gamma/self.beta)*self.ddp)
            #    for n in range(1, self.cL):
            #        drdt[n,:,:] = (
            #                lvln[n,:,:]
            #                - n*self.gamma*wig[n,:,:]
            #                + self.coef_up[n]*self.ddp[n+1,:,:]
            #                + self.coef_down[n]*n*theta_f[n-1,:,:]
            #                )
            if self.bath_type == "none":
                # Terms: Liouvillian
                drdt.insert(n,
                        lvln
                        )
            elif self.bath_type in {"debye","debye_correction3_cl"}:
                if self.cL ==0:
                    # Terms: Liouvillian + friction i.e. n = n
                    drdt.insert(n,
                            lvln
                            )
                else:
                    # Terms: Liouvillian + friction i.e. n = n
                    drdt.insert(n,
                            lvln
                            - np.dot(n,self.gammas)*wig
                            )
                    # rho_+k
                    if sum(n)>0:
                        for k in range(self.cK+1):
                            if n[k]!=0:
                                n_temp = list(n)
                                n_temp[k] -= 1
                                drdt.insert(n_temp,
                                        self.coef_up[n_temp[k],k]*self.ddp
                                        )
                    # rho_-k
                    if (sum(n)<self.cL) and (self.eta != 0):
                        # k > 0
                        for k in range(1,self.cK+1):
                            if n[k] != self.cL:
                                n_temp = list(n)
                                n_temp[k] += 1
                                drdt.insert(n_temp,
                                        self.coef_down[n_temp[k],k]
                                        * n_temp[k]*self.cs[k]*self.ddp
                                        )
                        # k = 0
                        if n[0] != self.cL:
                            n_temp = list(n)
                            n_temp[0] += 1
                            if self.bath_type=="debye":
                                drdt.insert(n_temp,
                                        self.coef_down[n_temp[0],0] * n_temp[0]
                                        * (self.eta*self.gamma/self.beta)
                                        *(self.ddp-self.beta*self.gamma
                                                            *self.q_mat*wig)
                                        )
                            elif self.bath_type=="debye_correction3_cl":
                                drdt.insert(n_temp,
                                        self.coef_down[n_temp[0],0] * n_temp[0]
                                        * (
                                            (self.eta*self.gamma/self.mass)
                                                                *self.p_mat*wig
                                            + (self.eta*self.gamma/self.beta)
                                                                *self.ddp
                                            )
                                        )
            else:
                raise SystemExit("Chosen bath type not supported for this sim")

    def write_density(self, filename, time):
        with open(filename, "a") as f:
            f.write(formatflt(time))
            data = self.dp*np.sum(self.ados[self.n0],axis=1)
            for i in range(self.n_q):
                f.write(formatflt(data[i]))
            f.write("\n")

    def check_stability(self, maxval):
        outcome = True
        if np.amax(self.dp*np.sum(self.ados[self.n0],axis=1))>maxval:
                outcome = False
        return outcome

    def write_ado_trend(self, filename, time):
        outcome = np.zeros([self.cL+1])
        for n, rho in self.ados.items():
            outcome[sum(n)] += np.amax(abs(self.ados[self.n0]))
        with open(filename, "a") as f:
            f.write(formatflt(time))
            for i in range(self.cL+1):
                f.write(formatflt(outcome[i]))
            f.write("\n")

    def postmultiply_by_qs(self):
        #self.ados[self.n0] *= self.q_mat
        for key in self.ados:
            self.ados[key] *= self.q_mat

    def write_tcf_qq(self, filename_real, filename_imag, time):
        write_complex_pair(filename_real, filename_imag, time,
            self.dq*self.dp*np.sum(self.ados[self.n0]*self.q_mat))

    def write_tcf_q2q2(self, filename_real, filename_imag, time):
        write_complex_pair(filename_real, filename_imag, time,
            self.dq*self.dp*np.sum(self.ados[self.n0]*self.q_mat*self.q_mat))

    def write_tcf_q4q4(self, filename_real, filename_imag, time):
        write_complex_pair(filename_real, filename_imag, time,
            self.dq*self.dp*np.sum(self.ados[self.n0]
                *self.q_mat*self.q_mat*self.q_mat*self.q_mat)
            )

 
#******************************************************************************
#******************************************************************************

class SimQHEOM_Mats_EE(SimBase,SimMats):
    """ Quantum HEOM simulation with Matsubara terms"""
    def __init__(self, inp, pot, state):
        SimBase.__init__(self,inp)
        SimMats.__init__(self,inp)
        self.n_EE = inp.n_EE
        if self.n_EE>len(self.qs):
            raise SystemExit("Cannot have more energy eigenstates than"
                    " point on the position grid")

        self.rho0 = np.zeros([self.n_EE,self.n_EE],dtype=complex)
        # Setting the initial state
        self.rho0[:,:],self.evals,self.evecs = state.rho_EE(self.qs,self.n_EE)
        # Normalising the initial state
        print("Normalisation constant")
        print(abs(np.trace(self.rho0)))
        self.rho0 /= np.trace(self.rho0)

        # Making Hamiltonian matrix
        self.ham = np.zeros([self.n_EE, self.n_EE])
        for i in range(self.n_EE):
            #print(self.evals[i])
            self.ham[i,i] = self.evals[i]
        self.iham = (-1j/hbar)*self.ham

        # Making q operator matrix
        self.q_mat = np.zeros([self.n_EE, self.n_EE])
        for i in range(self.n_EE):
            for j in range(self.n_EE):
                self.q_mat[i,j] = np.vdot(self.evecs[:,i], self.evecs[:,j]*self.qs)
        self.qq_mat = self.q_mat@self.q_mat
        self.iq_mat = (-1j/hbar)*self.q_mat

        self.dq_mat = self.evecs[:,:self.n_EE].T@dvr.derivative1(self.qs)@self.evecs[:,:self.n_EE]

        # Setting cutoff
        self.cutoff = np.amax(np.abs(self.rho0))*self.cutoff_coef

        # Work variables
        self.ados = Tensor()
        self.ados.insert(self.n0,self.rho0)
        self.k1 = Tensor()
        self.k2 = Tensor()
        self.k3 = Tensor()
        self.k4 = Tensor()
        self.rho_q = np.zeros([self.n_EE, self.n_EE], dtype=complex)
        self.q_rho = np.zeros([self.n_EE, self.n_EE], dtype=complex)

    def update_pot(self, pot):
        raise SystemExit("Updating potential is not supported for eigenvector "
                         "representation")

    def get_kubo_q_state(self, state):
        new_rho,new_evals,new_evecs = state.rho_EE_kubo_q(
                                                self.qs,self.n_EE,self.q_mat)
        self.ados = Tensor()
        self.ados.insert(self.n0,new_rho)
        self.k1 = Tensor()
        self.k2 = Tensor()
        self.k3 = Tensor()
        self.k4 = Tensor()

    def get_kubo_q2_state(self, state):
        new_rho,new_evals,new_evecs = state.rho_EE_kubo_q2(
                                                self.qs,self.n_EE,self.q_mat)
        self.ados = Tensor()
        self.ados.insert(self.n0,new_rho)
        self.k1 = Tensor()
        self.k2 = Tensor()
        self.k3 = Tensor()
        self.k4 = Tensor()

    def derivative(self, ados, coef, drdt):
        drdt.clear()
        for n, rho_raw in ados.items():
            rho = coef*rho_raw
            if self.bath_type == "none":
                # Terms: Liouvillian + friction i.e. n = n
                drdt.insert(n,
                            self.iham@rho - rho@self.iham
                        )
                # NB: iham includes -i/hbar
            elif self.bath_type in {"debye","debye_correction","debye_correction2","white"}:
                if self.bath_type=="debye":
                    # Terms: Liouvillian + friction i.e. n = n
                    drdt.insert(n,
                                self.iham@rho - rho@self.iham
                                - np.dot(n,self.gammas)*rho
                            )
                    # NB: iham includes -i/hbar
                elif self.bath_type in {"debye_correction"}:
                    drdt.insert(n,
                                self.iham@rho - rho@self.iham
                                - np.dot(n,self.gammas)*rho
                                - self.low_T_coef_R*(
                                                self.qq_mat@rho
                                                - 2*self.q_mat@rho@self.q_mat
                                                + rho@self.qq_mat)
                            )
                elif self.bath_type=="debye_correction2":
                    drdt.insert(n,
                                self.iham@rho - rho@self.iham
                                - np.dot(n,self.gammas)*rho
                                - (self.low_T_coef_R+self.low_T_coef_I)*self.qq_mat@rho
                                - (self.low_T_coef_R-self.low_T_coef_I)*rho@self.qq_mat
                                + 2*self.low_T_coef_R*self.q_mat@rho@self.q_mat
                            )
                elif self.bath_type in {"white"}:
                    # Note that the low_T_coef variables are used for coefficients
                    # of delta function terms in the bath TCF decomposition
                    drdt.insert(n,
                                self.iham@rho - rho@self.iham
                                - np.dot(n,self.gammas)*rho
                                - self.low_T_coef_R*(
                                                self.qq_mat@rho
                                                - 2*self.q_mat@rho@self.q_mat
                                                + rho@self.qq_mat
                                                )
                                - self.low_T_coef_I*(
                                                self.q_mat@self.dq_mat@rho
                                                +self.q_mat@rho@self.dq_mat
                                                -self.dq_mat@rho@self.q_mat
                                                -rho@self.dq_mat@self.q_mat
                                                #self.q_mat@self.dq_mat@rho
                                                #-self.q_mat@rho@self.dq_mat
                                                #-self.dq_mat@rho@self.q_mat
                                                #+rho@self.dq_mat@self.q_mat
                                    )
                            )

                # Calculating the "Commutator and anticommutator"
                # NB: iqs = -i/hbar * qs
                self.q_rho = self.iq_mat@rho
                self.rho_q = rho@self.iq_mat
                comm = self.q_rho-self.rho_q
                # rho_+k
                if sum(n)>0:
                    for k in range(self.cK+1):
                        if n[k]!=0:
                            n_temp = list(n)
                            n_temp[k] -= 1
                            drdt.insert(n_temp,self.coef_up[n_temp[k],k]*comm)
                # rho_-k
                if (sum(n)<self.cL) and (self.eta != 0):
                    # k > 0
                    for k in range(1,self.cK+1):
                        if n[k] != self.cL:
                            n_temp = list(n)
                            n_temp[k] += 1
                            drdt.insert(n_temp,
                                    self.coef_down[n_temp[k],k]
                                    * n_temp[k]*self.cs[k]*comm
                                    )
                    # k = 0
                    if n[0] != self.cL:
                        n_temp = list(n)
                        n_temp[0] += 1
                        drdt.insert(n_temp,
                                (self.coef_down[n_temp[0],0] * n_temp[0])
                                * (self.cs[0]*self.q_rho
                                    - self.cs[0].conjugate()*self.rho_q)
                                )
            else:
                raise SystemExit("Unknown bath_type selected")

    def write_density_EE(self, filename, time):
        with open(filename, "a") as f:
            f.write(formatflt(time))
            for i in range(self.n_EE):
                f.write(formatflt(np.real(self.ados[self.n0][i,i])))
            f.write("\n")

    def write_eval_grid(self,filename):
        with open(filename, "w+") as f:
            f.write(formatflt(0.0))
            for i in range(self.n_EE):
                f.write(formatflt(self.evals[i]))
            f.write("\n")

    def write_density(self, filename, time):
        # Constructing real space density
        diff = self.n_q-self.n_EE
        # Padding the density matrix with zeros
        density = np.pad(self.ados[self.n0], ((0,diff),(0,diff)),
                                                'constant', constant_values=0)
        # Converting the density matrix to position representation
        density = self.evecs@density@self.evecs.T.conj()
        with open(filename, "a") as f:
            f.write(formatflt(time))
            for i in range(self.n_q):
                f.write(formatflt(np.real(density[i,i])/self.dq))
            f.write("\n")

    def postmultiply_by_qs(self):
        #self.ados[self.n0] = self.ados[self.n0]@self.q_mat
        ## Multiplying all ados => gives nonsense
        #for n,rho in self.ados.items():
        #    rho = rho@self.q_mat
        #
        # It is true that all ADOs must be multiplied by q_mat
        for key in self.ados:
            self.ados[key]  = self.ados[key]@self.q_mat


    def write_tcf_qq(self, filename_real, filename_imag, time):
        write_complex_pair(filename_real, filename_imag, time,
            np.trace(self.ados[self.n0]@self.q_mat))

    def write_tcf_q2q2(self, filename_real, filename_imag, time):
        write_complex_pair(filename_real, filename_imag, time,
            np.trace(self.ados[self.n0]@self.q_mat@self.q_mat))

    def write_tcf_q4q4(self, filename_real, filename_imag, time):
        write_complex_pair(filename_real, filename_imag, time,
            np.trace(self.ados[self.n0]@self.q_mat@self.q_mat@self.q_mat@self.q_mat)
            )

    def check_stability(self, maxval):
        outcome = True
        for i in range(self.n_EE):
            if abs(self.ados[self.n0][i,i])>maxval:
                outcome = False
                break
        return outcome

    def write_ado_trend(self, filename, time):
        outcome = np.zeros([self.cL+1])
        for n, rho in self.ados.items():
            outcome[sum(n)] += np.amax(abs(self.ados[self.n0]))
        with open(filename, "a") as f:
            f.write(formatflt(time))
            for i in range(self.cL+1):
                f.write(formatflt(outcome[i]))
            f.write("\n")

#********************************************************************************
#********************************************************************************

class SimQHEOM_Mats_EE_heomc(SimBase,SimMats):
    """ Quantum HEOM simulation with Matsubara terms"""
    def __init__(self, inp, pot, state):
        SimBase.__init__(self,inp)
        SimMats.__init__(self,inp)
        self.n_EE = inp.n_EE
        if self.n_EE>len(self.qs):
            raise SystemExit("Cannot have more energy eigenstates than"
                    " point on the position grid")

        self.rho0 = np.zeros([self.n_EE,self.n_EE],dtype=complex)
        # Setting the initial state
        self.rho0[:,:],self.evals,self.evecs = state.rho_EE(self.qs,self.n_EE)
        # Normalising the initial state
        print("Normalisation constant")
        print(abs(np.trace(self.rho0)))
        self.rho0 /= np.trace(self.rho0)

        # Making Hamiltonian matrix
        self.ham = np.zeros([self.n_EE, self.n_EE])
        for i in range(self.n_EE):
            #print(self.evals[i])
            self.ham[i,i] = self.evals[i]
        self.iham = (-1j/hbar)*self.ham

        # Making q operator matrix
        self.q_mat = np.zeros([self.n_EE, self.n_EE])
        for i in range(self.n_EE):
            for j in range(self.n_EE):
                self.q_mat[i,j] = np.vdot(self.evecs[:,i], self.evecs[:,j]*self.qs)
        self.qq_mat = self.q_mat@self.q_mat
        self.iq_mat = (-1j/hbar)*self.q_mat

        self.dq_mat = self.evecs[:,:self.n_EE].T@dvr.derivative1(self.qs)@self.evecs[:,:self.n_EE]

        flattening = get_flattening_coefficients(self.cK, self.cL)

        #if "kubo" in inp.heom_switch:
        #    if inp.heom_switch=="kubo_qq":
        #        new_rho,new_evals,new_evecs = state.rho_EE_kubo_q(
        #                                        self.qs,self.n_EE,self.q_mat)
        #    elif inp.heom_switch=="kubo_q2q2":
        #        new_rho,new_evals,new_evecs = state.rho_EE_kubo_q2(
        #                                        self.qs,self.n_EE,self.q_mat)
        #    self.rho0 = new_rho
                
        low_T_coefs = np.zeros([2,1],dtype=complex)
        low_T_coefs[0] = self.low_T_coef_R
        low_T_coefs[1] = self.low_T_coef_I

        write_real_mat("tmp_q_mat",self.q_mat)
        write_real_mat("tmp_dq_mat",self.dq_mat)
        write_cplx_mat("tmp_ham",self.ham)
        if inp.heom_switch=="kubo_qq":
            output_rho,new_evals,new_evecs = state.rho_EE_kubo_q(
                                                self.qs,self.n_EE,self.q_mat)
        elif inp.heom_switch=="kubo_q2q2":
            output_rho,new_evals,new_evecs = state.rho_EE_kubo_q2(
                                                self.qs,self.n_EE,self.q_mat)
        else:
            output_rho = self.rho0
        write_cplx_mat("tmp_rho0",output_rho)
        write_real_mat("tmp_gammas",np.array(self.gammas))
        write_cplx_mat("tmp_low_T_coefs",low_T_coefs)
        write_cplx_mat("tmp_cs",np.array(self.cs))
        write_int_mat("tmp_flattening",flattening)

        # Setting cutoff
        self.cutoff = np.amax(np.abs(self.rho0))*self.cutoff_coef

        # Work variables
        self.ados = Tensor()
        self.ados.insert(self.n0,self.rho0)
        self.k1 = Tensor()
        self.k2 = Tensor()
        self.k3 = Tensor()
        self.k4 = Tensor()
        self.rho_q = np.zeros([self.n_EE, self.n_EE], dtype=complex)
        self.q_rho = np.zeros([self.n_EE, self.n_EE], dtype=complex)


    def update_pot(self, pot):
        raise SystemExit("Updating potential is not supported for eigenvector "
                         "representation")

    def get_kubo_q_state(self, state):
        new_rho,new_evals,new_evecs = state.rho_EE_kubo_q(
                                                self.qs,self.n_EE,self.q_mat)
        self.ados = Tensor()
        self.ados.insert(self.n0,new_rho)
        self.k1 = Tensor()
        self.k2 = Tensor()
        self.k3 = Tensor()
        self.k4 = Tensor()

    def get_kubo_q2_state(self, state):
        new_rho,new_evals,new_evecs = state.rho_EE_kubo_q2(
                                                self.qs,self.n_EE,self.q_mat)
        self.ados = Tensor()
        self.ados.insert(self.n0,new_rho)
        self.k1 = Tensor()
        self.k2 = Tensor()
        self.k3 = Tensor()
        self.k4 = Tensor()

    def derivative(self, ados, coef, drdt):
        drdt.clear()
        for n, rho_raw in ados.items():
            rho = coef*rho_raw
            if self.bath_type == "none":
                # Terms: Liouvillian + friction i.e. n = n
                drdt.insert(n,
                            self.iham@rho - rho@self.iham
                        )
                # NB: iham includes -i/hbar
            elif self.bath_type in {"debye","debye_correction","debye_correction2","white"}:
                if self.bath_type=="debye":
                    # Terms: Liouvillian + friction i.e. n = n
                    drdt.insert(n,
                                self.iham@rho - rho@self.iham
                                - np.dot(n,self.gammas)*rho
                            )
                    # NB: iham includes -i/hbar
                elif self.bath_type in {"debye_correction"}:
                    drdt.insert(n,
                                self.iham@rho - rho@self.iham
                                - np.dot(n,self.gammas)*rho
                                - self.low_T_coef_R*(
                                                self.qq_mat@rho
                                                - 2*self.q_mat@rho@self.q_mat
                                                + rho@self.qq_mat)
                            )
                elif self.bath_type=="debye_correction2":
                    drdt.insert(n,
                                self.iham@rho - rho@self.iham
                                - np.dot(n,self.gammas)*rho
                                - (self.low_T_coef_R+self.low_T_coef_I)*self.qq_mat@rho
                                - (self.low_T_coef_R-self.low_T_coef_I)*rho@self.qq_mat
                                + 2*self.low_T_coef_R*self.q_mat@rho@self.q_mat
                            )
                elif self.bath_type in {"white"}:
                    # Note that the low_T_coef variables are used for coefficients
                    # of delta function terms in the bath TCF decomposition
                    drdt.insert(n,
                                self.iham@rho - rho@self.iham
                                - np.dot(n,self.gammas)*rho
                                - self.low_T_coef_R*(
                                                self.qq_mat@rho
                                                - 2*self.q_mat@rho@self.q_mat
                                                + rho@self.qq_mat
                                                )
                                - self.low_T_coef_I*(
                                                self.q_mat@self.dq_mat@rho
                                                +self.q_mat@rho@self.dq_mat
                                                -self.dq_mat@rho@self.q_mat
                                                -rho@self.dq_mat@self.q_mat
                                                #self.q_mat@self.dq_mat@rho
                                                #-self.q_mat@rho@self.dq_mat
                                                #-self.dq_mat@rho@self.q_mat
                                                #+rho@self.dq_mat@self.q_mat
                                    )
                            )

                # Calculating the "Commutator and anticommutator"
                # NB: iqs = -i/hbar * qs
                self.q_rho = self.iq_mat@rho
                self.rho_q = rho@self.iq_mat
                comm = self.q_rho-self.rho_q
                # rho_+k
                if sum(n)>0:
                    for k in range(self.cK+1):
                        if n[k]!=0:
                            n_temp = list(n)
                            n_temp[k] -= 1
                            drdt.insert(n_temp,self.coef_up[n_temp[k],k]*comm)
                # rho_-k
                if (sum(n)<self.cL) and (self.eta != 0):
                    # k > 0
                    for k in range(1,self.cK+1):
                        if n[k] != self.cL:
                            n_temp = list(n)
                            n_temp[k] += 1
                            drdt.insert(n_temp,
                                    self.coef_down[n_temp[k],k]
                                    * n_temp[k]*self.cs[k]*comm
                                    )
                    # k = 0
                    if n[0] != self.cL:
                        n_temp = list(n)
                        n_temp[0] += 1
                        drdt.insert(n_temp,
                                (self.coef_down[n_temp[0],0] * n_temp[0])
                                * (self.cs[0]*self.q_rho
                                    - self.cs[0].conjugate()*self.rho_q)
                                )
            else:
                raise SystemExit("Unknown bath_type selected")

    def write_density_EE(self, filename, time):
        with open(filename, "a") as f:
            f.write(formatflt(time))
            for i in range(self.n_EE):
                f.write(formatflt(np.real(self.ados[self.n0][i,i])))
            f.write("\n")

    def write_eval_grid(self,filename):
        with open(filename, "w+") as f:
            f.write(formatflt(0.0))
            for i in range(self.n_EE):
                f.write(formatflt(self.evals[i]))
            f.write("\n")

    def write_density(self, filename, time):
        # Constructing real space density
        diff = self.n_q-self.n_EE
        # Padding the density matrix with zeros
        density = np.pad(self.ados[self.n0], ((0,diff),(0,diff)),
                                                'constant', constant_values=0)
        # Converting the density matrix to position representation
        density = self.evecs@density@self.evecs.T.conj()
        with open(filename, "a") as f:
            f.write(formatflt(time))
            for i in range(self.n_q):
                f.write(formatflt(np.real(density[i,i])/self.dq))
            f.write("\n")

    def postmultiply_by_qs(self):
        #self.ados[self.n0] = self.ados[self.n0]@self.q_mat
        ## Multiplying all ados => gives nonsense
        #for n,rho in self.ados.items():
        #    rho = rho@self.q_mat
        #
        # It is true that all ADOs must be multiplied by q_mat
        for key in self.ados:
            self.ados[key]  = self.ados[key]@self.q_mat


    def write_tcf_qq(self, filename_real, filename_imag, time):
        write_complex_pair(filename_real, filename_imag, time,
            np.trace(self.ados[self.n0]@self.q_mat))

    def write_tcf_q2q2(self, filename_real, filename_imag, time):
        write_complex_pair(filename_real, filename_imag, time,
            np.trace(self.ados[self.n0]@self.q_mat@self.q_mat))

    def write_tcf_q4q4(self, filename_real, filename_imag, time):
        write_complex_pair(filename_real, filename_imag, time,
            np.trace(self.ados[self.n0]@self.q_mat@self.q_mat@self.q_mat@self.q_mat)
            )

    def check_stability(self, maxval):
        outcome = True
        for i in range(self.n_EE):
            if abs(self.ados[self.n0][i,i])>maxval:
                outcome = False
                break
        return outcome

    def write_ado_trend(self, filename, time):
        outcome = np.zeros([self.cL+1])
        for n, rho in self.ados.items():
            outcome[sum(n)] += np.amax(abs(self.ados[self.n0]))
        with open(filename, "a") as f:
            f.write(formatflt(time))
            for i in range(self.cL+1):
                f.write(formatflt(outcome[i]))
            f.write("\n")

#********************************************************************************
#********************************************************************************
#class SimQHEOM_pyheom(SimBase,SimMats):
#    """ Quantum HEOM simulation with Matsubara terms"""
#    def __init__(self, inp, pot, state):
#        SimBase.__init__(self,inp)
#        SimMats.__init__(self,inp)
#        self.n_EE = inp.n_EE
#        if self.n_EE>len(self.qs):
#            raise SystemExit("Cannot have more energy eigenstates than"
#                    " points on the position grid")
#
#        self.rho0 = np.zeros([self.n_EE,self.n_EE],dtype=complex)
#        # Setting the initial state
#        self.rho0[:,:],self.evals,self.evecs = state.rho_EE(self.qs,self.n_EE)
#        # Normalising the initial state
#        print("Normalisation constant")
#        print(abs(np.trace(self.rho0)))
#        self.rho0 /= np.trace(self.rho0)
#
#        # Making Hamiltonian matrix
#        self.ham = np.zeros([self.n_EE, self.n_EE])
#        for i in range(self.n_EE):
#            #print(self.evals[i])
#            self.ham[i,i] = self.evals[i]
#
#        # Making q operator matrix
#        self.q_mat = np.zeros([self.n_EE, self.n_EE])
#        for i in range(self.n_EE):
#            for j in range(self.n_EE):
#                self.q_mat[i,j] = np.vdot(self.evecs[:,i], self.evecs[:,j]*self.qs)
#
#        self.J = pyheom.Drudian(self.eta,self.gamma)
#        if self.sop_decomposition_type=="matsubara":
#            self.corr_dict = pyheom.noise_decomposition(self.J, T=1/self.beta, type_LTC = 'MSD', n_MSD = self.cK)
#        elif self.sop_decomposition_type=="pade_N-1/N":
#            self.corr_dict = pyheom.noise_decomposition( self.J, T=1/self.beta, type_LTC = 'PSD', n_PSD = self.cK, type_PSD = "N-1/N")
#        elif self.sop_decomposition_type=="pade_N+1/N":
#            self.corr_dict = pyheom.noise_decomposition( self.J, T=1/self.beta, type_LTC = 'PSD', n_PSD = self.cK, type_PSD = "N+1/N")
#        elif self.sop_decomposition_type=="pade_N/N":
#            self.corr_dict = pyheom.noise_decomposition( self.J, T=1/self.beta, type_LTC = 'PSD', n_PSD = self.cK, type_PSD = "N/N")
#        else:
#            raise SystemExit("Unknown SOP decomposition scheme")
#        self.noises = [dict(V=self.q_mat, C=self.corr_dict)]
#        print("Creating HEOM object")
#        self.h = pyheom.HEOM(self.ham,
#                            self.noises,
#                            max_tier = self.cL,
#                            matrix_type='dense',
#                            #matrix_type='sparse',
#                            hierarchy_connection='loop',)
#        print("setting rho0")
#        self.h.set_rho(self.rho0)
#
##    def update_pot(self, pot):
##        raise SystemExit("Updating potential is not supported for eigenvector "
##                         "representation")
##
#    def get_kubo_q_state(self, state):
#        new_rho,new_evals,new_evecs = state.rho_EE_kubo_q(
#                                                self.qs,self.n_EE,self.q_mat)
#        self.h.set_rho(new_rho)
#
#    def get_kubo_q2_state(self, state):
#        new_rho,new_evals,new_evecs = state.rho_EE_kubo_q2(
#                                                self.qs,self.n_EE,self.q_mat)
#        self.h.set_rho(new_rho)
##
##    def derivative(self, ados, coef, drdt):
##        drdt.clear()
##        for n, rho_raw in ados.items():
##            rho = coef*rho_raw
##            if self.bath_type == 0:
##                # Terms: Liouvillian + friction i.e. n = n
##                drdt.insert(n,
##                            self.iham@rho - rho@self.iham
##                        )
##                # NB: iham includes -i/hbar
##            elif self.bath_type in {1,2,3}:
##                if self.bath_type==1:
##                    # Terms: Liouvillian + friction i.e. n = n
##                    drdt.insert(n,
##                                self.iham@rho - rho@self.iham
##                                - np.dot(n,self.gammas)*rho
##                            )
##                    # NB: iham includes -i/hbar
##                elif self.bath_type==2:
##                    drdt.insert(n,
##                                self.iham@rho - rho@self.iham
##                                - np.dot(n,self.gammas)*rho
##                                - self.low_T_coef_R*(
##                                                self.qq_mat@rho
##                                                - 2*self.q_mat@rho@self.q_mat
##                                                + rho@self.qq_mat)
##                            )
##                elif self.bath_type==3:
##                    drdt.insert(n,
##                                self.iham@rho - rho@self.iham
##                                - np.dot(n,self.gammas)*rho
##                                - (self.low_T_coef_R+self.low_T_coef_I)*self.qq_mat@rho
##                                - (self.low_T_coef_R-self.low_T_coef_I)*rho@self.qq_mat
##                                + 2*self.low_T_coef_R*self.q_mat@rho@self.q_mat
##                            )
##
##                # Calculating the "Commutator and anticommutator"
##                # NB: iqs = -i/hbar * qs
##                self.q_rho = self.iq_mat@rho
##                self.rho_q = rho@self.iq_mat
##                comm = self.q_rho-self.rho_q
##                # rho_+k
##                if sum(n)>0:
##                    for k in range(self.cK+1):
##                        if n[k]!=0:
##                            n_temp = list(n)
##                            n_temp[k] -= 1
##                            drdt.insert(n_temp,self.coef_up[n_temp[k],k]*comm)
##                # rho_-k
##                if (sum(n)<self.cL) and (self.eta != 0):
##                    # k > 0
##                    for k in range(1,self.cK+1):
##                        if n[k] != self.cL:
##                            n_temp = list(n)
##                            n_temp[k] += 1
##                            drdt.insert(n_temp,
##                                    self.coef_down[n_temp[k],k]
##                                    * n_temp[k]*self.cs[k]*comm
##                                    )
##                    # k = 0
##                    if n[0] != self.cL:
##                        n_temp = list(n)
##                        n_temp[0] += 1
##                        drdt.insert(n_temp,
##                                (self.coef_down[n_temp[0],0] * n_temp[0])
##                                * (self.cs[0]*self.q_rho
##                                    - self.cs[0].conjugate()*self.rho_q)
##                                )
##            else:
##                raise SystemExit("Unknown bath_type selected")
##
#    def write_density_EE(self, filename, time, rho):
#        with open(filename, "a") as f:
#            f.write(formatflt(time))
#            for i in range(self.n_EE):
#                f.write(formatflt(np.real(rho[i,i])))
#            f.write("\n")
#
#    def write_eval_grid(self,filename):
#        with open(filename, "w+") as f:
#            f.write(formatflt(0.0))
#            for i in range(self.n_EE):
#                f.write(formatflt(self.evals[i]))
#            f.write("\n")
#
#    def write_density(self, filename, time, rho):
#        # Constructing real space density
#        diff = self.n_q-self.n_EE
#        # Padding the density matrix with zeros
#        density = np.pad(rho, ((0,diff),(0,diff)),
#                                                'constant', constant_values=0)
#        # Converting the density matrix to position representation
#        density = self.evecs@density@self.evecs.T.conj()
#        with open(filename, "a") as f:
#            f.write(formatflt(time))
#            for i in range(self.n_q):
#                f.write(formatflt(np.real(density[i,i])/self.dq))
#            f.write("\n")
#
#    def postmultiply_by_qs(self):
#        #self.ados[self.n0] = self.ados[self.n0]@self.q_mat
#        ## Multiplying all ados => gives nonsense
#        #for n,rho in self.ados.items():
#        #    rho = rho@self.q_mat
#        #
#        # It is true that all ADOs must be multiplied by q_mat
#        for i in range(self.h.rho_h.shape[2]):
#            self.h.rho_h[:,:,i] = self.h.rho_h[:,:,i]@self.q_mat
#            #self.ados[key]  = self.ados[key]@self.q_mat
#
#
#    def write_tcf_qq(self, filename_real, filename_imag, time, rho):
#        write_complex_pair(filename_real, filename_imag, time,
#            np.trace(rho@self.q_mat))
#
#    def write_tcf_q2q2(self, filename_real, filename_imag, time, rho):
#        print("Before I even start doing stuff")
#        write_complex_pair(filename_real, filename_imag, time,
#            np.trace(rho@self.q_mat@self.q_mat))
#
#    def write_tcf_q4q4(self, filename_real, filename_imag, time, rho):
#        write_complex_pair(filename_real, filename_imag, time,
#            np.trace(rho@self.q_mat@self.q_mat@self.q_mat@self.q_mat))
##
##    def check_stability(self, maxval):
##        outcome = True
##        for i in range(self.n_EE):
##            if abs(self.ados[self.n0][i,i])>maxval:
##                outcome = False
##                break
##        return outcome
##
##    def write_ado_trend(self, filename, time):
##        outcome = np.zeros([self.cL+1])
##        for n, rho in self.ados.items():
##            outcome[sum(n)] += np.amax(abs(self.ados[self.n0]))
##        with open(filename, "a") as f:
##            f.write(formatflt(time))
##            for i in range(self.cL+1):
##                f.write(formatflt(outcome[i]))
##            f.write("\n")
##
