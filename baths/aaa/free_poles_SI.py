#!/usr/bin/env python3
import numpy as np
import sys
import matplotlib.pyplot as plt
from scipy.linalg import eig

# This file contains the coefficients for the poles of the bath correlation function

data=[[1,-0.007050248,0.002482056,188.7752364,-164.9169644],
[2,0.102762949,0.116652539,114.5230822,-130.046728],
[3,0.581303135,-0.644477046,69.37948966,-100.2912245],
[4,-1.455575575,-1.650266045,40.55914002,-74.97480723],
[5,-2.988444895,1.187910091,22.52508418,-54.17566836],
[6,-0.622081237,3.115931332,11.81278425,-37.82551602],
[7,1.30388447,2.013121889,5.852688004,-25.57526889],
[8,1.382279145,0.625686611,2.752595646,-16.83290689],
[9,0.8611309,0.016022,1.241895418,-10.8518732],
[10,0.443366368,-0.118339778,0.548086947,-6.885813753],
[11,0.212002544,-0.101557778,0.243348806,-4.31568345],
[12,0.09893136,-0.063028266,0.113960488,-2.680510104],
[13,0.046158834,-0.034415289,0.059450913,-1.654822486],
[14,0.021773639,-0.017637193,0.033910866,-1.017022414],
[15,0.010348793,-0.008760985,0.019524006,-0.62200475],
[16,0.004904562,-0.004276533,0.010782038,-0.378232958],
[17,0.002297321,-0.002058266,0.005782001,-0.228876357],
[18,0.001065577,-0.00097944,0.003293429,-0.138297459],
[19,0.00049681,-0.000464088,0.002202629,-0.083659132],
[20,0.000235641,-0.000219394,0.001675626,-0.050533949],
[21,0.000112977,-0.000102543,0.001268055,-0.03029634],
[22,5.37E-05,-4.70E-05,0.000861682,-0.017943686],
[23,2.49E-05,-2.12E-05,0.000496313,-0.010482398],
[24,1.12E-05,-9.49E-06,0.000227853,-0.006041449],
[25,4.85E-06,-4.23E-06,6.86E-05,-0.003437852],
[26,2.05E-06,-1.89E-06,-9.00E-06,-0.001928643],
[27,8.36E-07,-8.31E-07,-3.82E-05,-0.001060056],
[28,3.24E-07,-3.56E-07,-4.17E-05,-0.000568788],
[29,1.21E-07,-1.49E-07,-3.41E-05,-0.000297603],
[30,4.48E-08,-6.28E-08,-2.43E-05,-0.000147216],
[31,1.62E-08,-2.71E-08,-1.54E-05,-5.76E-05]]
# Convert data to numpy array
data = np.array(data)
# Extract columns
K = data.shape[0]  # Number of poles
K = 29 # Number of poles
n = data[:K, 0].astype(int)  # Pole number
eta = data[:K, 1] + 1.j* data[:K,2]            # Coefficient residue
xi = data[:K, 3] + 1.j* data[:K,4]             # Pole frequency

# define S(omega) and the functions S_plus and S_minus
def S(omega):
    return np.sum(eta[:,np.newaxis] / (omega[np.newaxis,:] - xi[:,np.newaxis]), axis=0) + np.sum(eta[:,np.newaxis].conj() / (omega[np.newaxis,:] - xi[:,np.newaxis].conj()), axis=0)
def S_plus(omega): return S(omega) + S(-omega)
def S_minus(omega): return S(omega) - S(-omega)
def A(omega): return S_plus(omega) / S_minus(omega)

### Calculate the poles and residues of A(ω) = S⁺(ω)/S⁻(ω) using the generalized eigenvalue problem ###
# Make the x_ks, z_ks, f_ks arrays for the pole and residue calculation
x_ks = np.zeros(4*K, dtype=complex)
z_ks = np.zeros(4*K, dtype=complex)
f_ks = np.zeros(4*K, dtype=complex)
for ii in range(K):
    idx=4*ii
    x_ks[idx] = eta[ii]; x_ks[idx+1] = eta[ii]; x_ks[idx+2] = eta[ii].conj(); x_ks[idx+3] = eta[ii].conj()
    z_ks[idx] = xi[ii]; z_ks[idx+1] = -xi[ii]; z_ks[idx+2] = xi[ii].conj(); z_ks[idx+3] = -xi[ii].conj()
    f_ks[idx] = 1; f_ks[idx+1] = -1; f_ks[idx+2] = 1; f_ks[idx+3] = -1
# make the E and B matrices
E=np.zeros((4*K+1,4*K+1), dtype=complex)
B=np.eye(4*K+1, dtype=complex); B[0,0]=0
E[0,1:]=x_ks
E[1:,0]=1
E[1:,1:]=np.diag(z_ks)
# Solve the generalized eigenvalue problem
A_poles, _ = eig(E, B)
A_poles = A_poles[~np.isinf(A_poles)]  # remove infs if any
# round the poles to remove numerical noise
A_poles = np.round(A_poles, decimals=10)
# remove the non pure imaginary poles (numerical noise)
# A_poles = A_poles[np.abs(A_poles.real) < 1e-0]

# Calculate the residues using contour integration trick
A_res = np.zeros_like(A_poles, dtype=complex)  # residues of A(ω)
for k in range(len(A_poles)):
    dz = 1e-5 * np.exp(2j * np.pi * np.arange(40) / 40)  # small circle around each pole
    z_plus_dz = A_poles[k] + dz
    vals = A(z_plus_dz)
    A_res[k] = np.dot(vals, dz) / 40  # trapezoidal integration (approx.)
#remove residues that are too small (numerical noise)
A_res = A_res[np.abs(A_res) > 1e-5]
A_poles = A_poles[np.abs(A_res) > 1e-5]

# Define a function that recalculates A(ω) from the poles and residues
def A_from_poles(omega):
    return np.sum(A_res[:,np.newaxis] / (omega[np.newaxis,:]- A_poles[:,np.newaxis]), axis=0)


if(1):
    # Plot the fit of Aw from poles vs Aw from S(ω)
    fig, ax = plt.subplots(1, 1, figsize=(6.5, 3.2))  # width, height in inches
    # Set font size to match LaTeX \normalsize (about 10pt) 
    plt.rcParams.update({'font.size': 10})
    omega = np.linspace(-3000, 3000, 10000)  # Frequency range
    ax.plot(omega, A(omega).real, label=r'$A(\omega)$ from $S(\omega)$', color='blue', linestyle='-')
    ax.plot(omega, A_from_poles(omega).real, label=r'$A(\omega)$ from poles', color='red', linestyle='--')
    ax.plot(omega, A_from_poles(omega).imag, label=r'$A(\omega)$ from poles', color='green', linestyle='--')
    ax.set_xlabel('Frequency (ω)')
    ax.set_ylabel(r'$A(\omega)$')


# plot the poles of A(W)
if(1):
    fig, ax = plt.subplots(1, 1, figsize=(6.5, 3.2))  # width, height in inches
    plt.rcParams.update({'font.size': 10})
    ax.scatter(A_poles.real, A_poles.imag, color='blue', marker='o', label='Poles of $A(\omega)$')
    ax.axhline(0, color='black', linewidth=0.5, linestyle='--')
    ax.axvline(0, color='black', linewidth=0.5, linestyle='--')
    ax.set_xlabel('Real Part')
    ax.set_ylabel('Imaginary Part')
    ax.set_title('Poles of $A(\\omega)$ in the Complex Plane')
    # label each pole with its residue
    for i in range(len(A_poles)):
        ax.text(A_poles[i].real, A_poles[i].imag, f'{np.round(A_res[i],2):.2e}', fontsize=8, ha='right')
    ax.legend()


if(0):
    # plot S(w)
    fig, ax = plt.subplots(1, 1, figsize=(6.5, 3.2))  # width, height in inches
    omega = np.linspace(-300, 300, 10000)  # Frequency range
    ax.plot(omega, S(omega).real, label=r'$S(\omega)$', color='blue', linestyle='-')
    ax.plot(omega, S(omega).imag, label=r'$S(\omega)$', color='red', linestyle='--')
    ax.set_xlabel('Frequency (ω)')
    ax.set_ylabel(r'$S(\omega)$')


plt.show()
sys.exit()