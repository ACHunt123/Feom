#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

### Load the aaa files
folder='.files'
repoles = np.loadtxt(f'{folder}/pol_real.txt')
impoles = np.loadtxt(f'{folder}/pol_imag.txt')
poles = repoles + 1.j * impoles
reres = np.loadtxt(f'{folder}/res_real.txt')
imres = np.loadtxt(f'{folder}/res_imag.txt')
res = reres + 1.j * imres

# Setup the pole plot
fig, (pole_ax,func_ax) =  plt.subplots(1, 2, figsize=(12, 6))
pole_ax.set_title('Poles')
pole_ax.set_xlabel('Real Part')
pole_ax.set_ylabel('Imaginary Part')
pole_ax.grid(True)
pole_ax.plot(repoles, impoles, 'o', label=f'AAA, N={len(poles)}', color='blue')


# Frequency range
w = np.linspace(-100, 100, 2000)  # Adjust

# Calculate the radius of gyration exactly
beta=1 ; hbar=1 ; m=1 ; u = beta * hbar * w
def coth(x): return 1 / np.tanh(x)   
S_exact_values = (1/(m * beta * w**2)) * ((u/2) * coth(u/2) - 1)


# Define S(omega), which recreates the  function from the poles
def S_aaa(omega, res, poles):
    summ = 0+0j
    for i in range(len(poles)):
        summ += res[i] / (omega*(1+0j) - poles[i])
    return summ
S_aaa_values = np.array([S_aaa(omega, res, poles) for omega in w])  
# take out the constant term
c = np.average(S_aaa_values- S_exact_values)
S_aaa_values -= c

# Calculate the radius of gyration using N Matsubara frequencies/rp modes
N= 25 #ODD number of beads
mu = (N-1)/2 # Number of pairs of modes
ks = np.arange(1,  mu+1) # modes indices
wks = 2*N/(beta*hbar) *np.sin(ks * np.pi / (N))
pole_ax.plot(np.zeros_like(wks), wks, 'x', label=f'RP, N={N-1}', color='orange')
pole_ax.plot(np.zeros_like(wks), -wks, 'x', color='orange')

def S_rp(omega, wks):
    summ = 0+0j
    for wk in wks:
        summ += (2/(m*beta))*  1/(omega**2*(1+0j) + wk**2)
    return summ
S_rp_values = np.array([S_rp(omega, wks) for omega in w])

# Plot the results
func_ax.plot(w, S_aaa_values.real, label='Re[S_aaa(ω)]', color='red')
func_ax.plot(w, S_aaa_values.imag, label='Im[S_aaa(ω)]', color='green')
func_ax.plot(w, S_rp_values.real, label='Re[S_rp(ω)]', color='orange', linestyle='--')
func_ax.plot(w, S_rp_values.imag, label='Im[S_rp(ω)]', color='purple', linestyle='--')
# Plot the exact S(ω) for comparison

func_ax.plot(w, S_exact_values, label='Exact S(ω)', color='black', linestyle='--')
#plot the difference between the exact and the poles
func_ax.plot(w, S_aaa_values.real - S_exact_values, label='Re[S(ω)] - Exact', color='blue', linestyle=':')
func_ax.set_xlabel('Frequency (ω)')     
func_ax.set_ylabel('S(ω)')
func_ax.set_title('Bath Correlation Function S(ω)')
func_ax.legend()
pole_ax.legend()
# plt.ylim(-20, 20)
func_ax.set_xlim(-100, 100)
func_ax.grid()
plt.show()
