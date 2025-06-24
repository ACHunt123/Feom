#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
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
n = data[:, 0].astype(int)  # Pole number
eta = data[:, 1] + 1.j* data[:,2]               # Coefficient residue
xi = data[:, 3] + 1.j* data[:,4]             # Pole frequency

# plot the poles
plt.figure(figsize=(8, 6))
# plt.scatter(xi.real, xi.imag, color='blue', label='Poles')
# plt.scatter(xi.real, -xi.imag, color='blue')  # Mirror image for

# define S(omega)
def S(omega):
    one = np.ones_like(eta)
    return np.sum(eta / (omega*one - xi)) + np.sum(eta.conj() / (omega*one - xi.conj()))

w = np.linspace(-1000, 1000, 10000)  # Frequency range
S = np.array([S(omega) for omega in w])  # Calculate S(ω) for each frequency
plt.plot(w, S.real, label='Re[S(ω)]', color='red')
plt.plot(w, S+np.flip(S), label='S(ω) + S(-ω)', color='blue')
plt.plot(w, S-np.flip(S), label='S(ω) - S(-ω)', color='orange')
plt.plot(w, (S+np.flip(S))/(S-np.flip(S)), label='(S(ω) + S(-ω)) / (S(ω) - S(-ω))', color='purple')
# plt.plot(w, S.imag, label='Im[S(ω)]', color='green')
plt.axhline(0, color='black', lw=0.5, ls='--')
plt.axvline(0, color='black', lw=0.5, ls='--')
plt.xlabel('Frequency (ω)')
plt.ylabel('S(ω)')
plt.legend()
plt.ylim(-5, 5)
plt.show()