#!/usr/bin/env python
# File: general.py
""" General module with constants formatting etc."""
import numpy as np
import time

pi = 3.141592653589793238462643383279502884197169399375105
hbar = 1.000000000000000
au2energy = 4.359744650e-18
au2charge = 1.6021766208e-19 

def printobj(obj):
    for item in obj.__dict__:
        print(item, obj.__dict__[item])

def formatflt(flt):
    #return np.format_float_scientific(flt, precision=15, exp_digits=3)+"     "
    return "{: 20.16e}     ".format(flt)

def formatfltshort(flt):
    #return np.format_float_scientific(flt, precision=15, exp_digits=3)+"     "
    return "{: 6.4e}     ".format(flt)

def write_complex_pair(filename_real, filename_imag, time, value):
    with open(filename_real, "a") as re, open(filename_imag, "a") as im:
        re.write(formatflt(time))
        re.write(formatflt(np.real(value)))
        re.write("\n")
        im.write(formatflt(time))
        im.write(formatflt(np.imag(value)))
        im.write("\n")

class Timer(object):
    def __init__(self):
        self.time0 = time.time()
        self.process_time0 = time.process_time()
    def start(self):
        self.time0 = time.time()
        self.process_time0 = time.process_time()
    def end(self):
        self.time1 = time.time()
        self.process_time1 = time.process_time()
    def process(self):
        return self.process_time1-self.process_time0
    def time(self):
        return self.time1-self.time0
    def clock(self):
        return self.process_time1-self.process_time0

def write_int_mat(filename,mat):
    np.savetxt(filename,mat,fmt='%1d')
def write_real_mat(filename,mat):
    np.savetxt(filename,mat,fmt='%.18e')
def write_cplx_mat(filename,inpmat):
    mat = np.array(inpmat,dtype=complex)
    fmt_root = '(%.18e,%.18e)'
    if len(mat.shape)==1:
        n_cols = 1
    else:
        n_cols = mat.shape[1]
    fmt = fmt_root
    for i in range(n_cols-1):
        fmt += ' '+fmt_root
    np.savetxt(filename,mat,fmt=fmt)

def get_flattening_coefficients(K,L):
    from scipy.special import binom
    mat = np.zeros([L+1,K+1],dtype=int)
    for k in range(0,K+1):
        for n in range(0,L+1):
            mat[n,k]=binom(n+k,k+1)
    return mat


