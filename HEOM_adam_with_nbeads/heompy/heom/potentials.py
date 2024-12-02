#!/usr/bin/env python
# File: potentials.py
"""Potential object module."""
import math as m
import numpy as np

#******************************************************************************
#******************************************************************************

class PotEff(object):
    """Effective potential including renormalisation"""
    def __init__(self, pot, a_ren):
        self.pot = pot
        self.a_ren = a_ren

    def calc(self,x):
        return self.pot.calc(x) + 0.5*self.a_ren*x**2

    def switch(self):
        self.pot.switch()

    def diff(self,x):
        return self.pot.diff(x)+ self.a_ren*x

#******************************************************************************
#******************************************************************************

class PotHO(object):
    """Harmonic oscillator potential"""
    def __init__(self, inp):
        self.omega = inp.ho_omega
        self.mass = inp.mass
        self.k = self.mass*(self.omega**2)
        self.x0 = inp.ho_x0
        self.x0_switch = inp.ho_x0_switch
        self.const = inp.ho_const
        self.const_switch = inp.ho_const_switch

    def calc(self,x):
        return 0.5*self.k*(x-self.x0)**2 + self.const

    def switch(self):
        self.x0 = self.x0_switch
        self.const = self.const_switch

    def diff(self,x):
        return self.k*(x-self.x0)

#******************************************************************************
#******************************************************************************

class PotTani(object):
    """Potential from Tanimura(1991) article"""
    def __init__(self, inp):
        self.omega_a = inp.tani_omega_a
        self.f_coef = inp.tani_f_coef
        self.f_coef_switch = inp.tani_f_coef_switch
        self.mass = inp.mass

    def calc(self,x):
        return (0.5*self.mass*(self.omega_a**2)*(x**2)
                                        - self.f_coef*x)
    def switch(self):
        self.f_coef = self.f_coef_switch

    def diff(self,x):
        return self.mass*(self.omega_a**2)*x - self.f_coef

#******************************************************************************
#******************************************************************************

class PotMorse(object):
    """Morse potential"""
    def __init__(self, inp):
        self.alpha = inp.morse_alpha
        self.de = inp.morse_de
        self.x0 = inp.morse_x0
        self.x0_switch = inp.morse_x0_switch
        self.const = inp.morse_const
        self.const_switch = inp.morse_const_switch

    def calc(self,x):
        return ( self.de*( 1-np.exp(-self.alpha*(x-self.x0)) )**2
                + self.const)

    def switch(self):
        self.x0 = self.x0_switch
        self.const = self.const_switch

    def diff(self,x):
        return (2*self.de*self.alpha*( 1-np.exp(-self.alpha*(x-self.x0)) )
                * np.exp(self.alpha*(x-self.x0)) )

#******************************************************************************
#******************************************************************************

class PotPoly(object):
    """General polynomial potential"""
    def __init__(self, inp):
        self.coefs = inp.poly_coefs
        self.coefs_switch = inp.poly_coefs_switch
        self.x0 = inp.poly_x0
        self.x0_switch = inp.poly_x0_switch

    def calc(self,x):
        V = 0.0
        for i in range(len(self.coefs)):
            V += self.coefs[i]*(x-self.x0)**i
        return V

    def switch(self):
        self.x0 = self.x0_switch
        self.coefs = self.coefs_switch

    def diff(self,x):
        dV = 0.0
        for i in range(1,len(self.coefs)):
            dV += self.coefs[i]*(i)*(x-self.x0)**(i-1)
        return dV

#******************************************************************************
#******************************************************************************

class PotHydOH(object):
    """Raz's anharmonic potential inspired by hydrated OH from Joel Bowman's
        article"""
    def __init__(self, inp):
        self.c = inp.hydoh_c

    def calc(self,x):
        V = 0.1297*x**2 + self.c*(0.1657*x**4 - 0.2467*x**3)
        return V

    def switch(self):
        pass

    def diff(self,x):
        dV = 0.1297*2*x + self.c*(0.1657*4*x**3 - 0.2467*3*x**2)
        return dV 

#******************************************************************************
#******************************************************************************

class PotHarmonic(object):
    """Harmonic oscillator potential"""
    def __init__(self, inp):
        self.omega = inp.ho_omega
        self.mass = inp.mass
        self.k = self.mass*(self.omega**2)
        self.x0 = inp.ho_x0
        self.x0_switch = inp.ho_x0_switch
        self.const = inp.ho_const
        self.const_switch = inp.ho_const_switch

    def calc(self,x):
        return 0.5*self.k*(x-self.x0)**2 + self.const

    def switch(self):
        self.x0 = self.x0_switch
        self.const = self.const_switch

    def diff(self,x):
        return self.k*(x-self.x0)

#******************************************************************************
#******************************************************************************

class PotAnharmonic(object):
    """Harmonic oscillator potential"""
    def __init__(self, inp):
        self.omega = inp.ho_omega
        self.mass = inp.mass
        self.k = self.mass*(self.omega**2)
        self.x0 = inp.ho_x0
        self.x0_switch = inp.ho_x0_switch
        self.const = inp.ho_const
        self.const_switch = inp.ho_const_switch

    def calc(self,x):
        return self.k*(0.5*(x-self.x0)**2 + 0.1*(x-self.x0)**3
                                          + 0.01*(x-self.x0)**4) + self.const

    def switch(self):
        self.x0 = self.x0_switch
        self.const = self.const_switch

    def diff(self,x):
        return self.k*((x-self.x0) + 0.3*(x-self.x0)**2
                                          + 0.04*(x-self.x0)**3) 
#******************************************************************************
#******************************************************************************

class PotQuartic(object):
    """Harmonic oscillator potential"""
    def __init__(self, inp):
        self.omega = inp.ho_omega
        self.mass = inp.mass
        self.k = self.mass*(self.omega**2)
        self.x0 = inp.ho_x0
        self.x0_switch = inp.ho_x0_switch
        self.const = inp.ho_const
        self.const_switch = inp.ho_const_switch

    def calc(self,x):
        return 0.25*self.k*(x-self.x0)**4 + self.const

    def switch(self):
        self.x0 = self.x0_switch
        self.const = self.const_switch

    def diff(self,x):
        return self.k*(x-self.x0)**3

#******************************************************************************
#******************************************************************************

class PotHarmoquartic(object):
    """Harmonic oscillator potential"""
    def __init__(self, inp):
        self.omega = inp.ho_omega
        self.mass = inp.mass
        self.k = self.mass*(self.omega**2)
        self.x0 = inp.ho_x0
        self.x0_switch = inp.ho_x0_switch
        self.const = inp.ho_const
        self.const_switch = inp.ho_const_switch

    def calc(self,x):
        return self.k*(0.5*(x-self.x0)**2 + 5*(x-self.x0)**4) + self.const

    def switch(self):
        self.x0 = self.x0_switch
        self.const = self.const_switch

    def diff(self,x):
        return self.k*((x-self.x0) + 20*(x-self.x0)**3)

#******************************************************************************
#******************************************************************************

class PotChampagne(object):
    """1D OH bond-like potential"""
    def __init__(self, inp):
        self.mass = inp.mass
        self.d0 = 0.18748
        self.alpha = 1.1605
        self.x0 = inp.ho_x0
        self.x0_switch = inp.ho_x0_switch
        self.const = inp.ho_const
        self.const_switch = inp.ho_const_switch

    def calc(self,x):
        return self.d0*((1-np.exp(-self.alpha*(x-self.x0)))**2) + self.const

    def switch(self):
        self.x0 = self.x0_switch
        self.const = self.const_switch

    def diff(self,x):
        tmp = np.exp(-self.alpha*(x-self.x0))
        return 2*self.alpha*self.d0*tmp*(1-tmp)

#******************************************************************************
#******************************************************************************

class PotDW(object):
    def __init__(self, mass, eb, omega_b):
        self.mass = mass
        self.eb = eb
        self.omega_b = omega_b
        self.a1 = 0.5*mass*(omega_b**2)
        self.a2 = (mass**2)*(omega_b**4)/(16*eb)

    def calc(self,x):
        return -self.a1*(x**2) + self.a2*(x**4)

    def diff(self,x):
        return -2*self.a1*x + 4*self.a2*(x**3)

class PotDW1(PotDW):
    """DW1 potential from Topaler, Makri, J. Chem. Phys. 101, 7500 (1994)"""
    def __init__(self, inp):
        # energy barrier = 2085 cm^-1
        eb = 9.499959002321540e-3
        # barrier frequency omega_b = 500 cm^-1
        omega_b = 2.2781676264560e-3
        # originally the mass of the proton = 1836.15267344 a.u.
        mass = inp.mass
        PotDW.__init__(self,mass,eb,omega_b)

class PotDW2(PotDW):
    """DW2 potential from Topaler, Makri, J. Chem. Phys. 101, 7500 (1994)"""
    def __init__(self, inp):
        # energy barrier = 1043 cm^-1
        eb = 4.752257668787230e-3
        # barrier frequency omega_b = 500 cm^-1
        omega_b = 2.2781676264560e-3
        # originally the mass of the proton = 1836.15267344 a.u.
        mass = inp.mass
        PotDW.__init__(self,mass,eb,omega_b)
