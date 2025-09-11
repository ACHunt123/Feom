#!/usr/bin/env python3
import numpy as np
import sys
import matplotlib.pyplot as plt
from scipy.linalg import eig

# This file contains the coefficients for the poles of the bath correlation function
import numpy as np
from scipy.linalg import eig

class SpectralFunction:
    def __init__(self, poles=None, residues=None, K=None):
        if poles is None or residues is None:
            self._load_default_data(K)
        else:
            self.K = len(poles)
            self.eta = np.array(residues)
            self.xi = np.array(poles)
        
        self.calculate_A_poles_and_residues()

    def _load_default_data(self, K_override=None):
        # (Truncated default data for brevity — use full table in real code)
        data=np.array([[1,-0.007050248,0.002482056,188.7752364,-164.9169644],
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
        [31,1.62E-08,-2.71E-08,-1.54E-05,-5.76E-05]])
        K = K_override or data.shape[0]
        self.K = K
        self.eta = data[:K,1] + 1j*data[:K,2]
        self.xi = data[:K,3] + 1j*data[:K,4]

    def S(self, omega):
        term1 = self.eta[:, None] / (omega[None, :] - self.xi[:, None])
        term2 = self.eta[:, None].conj() / (omega[None, :] - self.xi[:, None].conj())
        return np.sum(term1 + term2, axis=0)

    def S_plus(self, omega):
        return self.S(omega) + self.S(-omega)

    def S_minus(self, omega):
        return self.S(omega) - self.S(-omega)

    def A(self, omega):
        return self.S_plus(omega) / self.S_minus(omega)

    def calculate_A_poles_and_residues(self):
        K = self.K
        eta = self.eta
        xi = self.xi
        x_ks = np.zeros(4*K, dtype=complex)
        z_ks = np.zeros(4*K, dtype=complex)
        f_ks = np.zeros(4*K, dtype=complex)

        for ii in range(K):
            idx = 4*ii
            x_ks[idx:idx+4] = [eta[ii], eta[ii], eta[ii].conj(), eta[ii].conj()]
            z_ks[idx:idx+4] = [xi[ii], -xi[ii], xi[ii].conj(), -xi[ii].conj()]
            f_ks[idx:idx+4] = [1, -1, 1, -1]

        E = np.zeros((4*K+1, 4*K+1), dtype=complex)
        B = np.eye(4*K+1, dtype=complex); B[0,0] = 0
        E[0, 1:] = x_ks
        E[1:, 0] = 1
        E[1:, 1:] = np.diag(z_ks)

        poles_raw, _ = eig(E, B)
        poles_clean = np.round(poles_raw[~np.isinf(poles_raw)], 10)

        # Store for later use
        self.A_poles = poles_clean
        self.A_res = self._calculate_residues(self.A_poles)

    def _calculate_residues(self, poles):
        dz = 1e-5 * np.exp(2j * np.pi * np.arange(40) / 40)
        residues = np.zeros_like(poles, dtype=complex)
        for k in range(len(poles)):
            z_plus_dz = poles[k] + dz
            vals = self.A(z_plus_dz)
            residues[k] = np.dot(vals, dz) / 40
        mask = np.abs(residues) > 1e-5
        self.A_poles = self.A_poles[mask]
        return residues[mask]

    def A_from_poles(self, omega):
        return np.sum(self.A_res[:, None] / (omega[None, :] - self.A_poles[:, None]), axis=0)

if __name__ == "__main__":
    sf = SpectralFunction()

    if(1):
        # Plot the fit of Aw from poles vs Aw from S(ω)
        fig, ax = plt.subplots(1, 1, figsize=(6.5, 3.2))  # width, height in inches
        # Set font size to match LaTeX \normalsize (about 10pt) 
        plt.rcParams.update({'font.size': 10})
        omega = np.linspace(-3000, 3000, 10000)  # Frequency range
        ax.plot(omega, sf.A(omega).real, label=r'$A(\omega)$ from $S(\omega)$', color='blue', linestyle='-')
        ax.plot(omega, sf.A_from_poles(omega).real, label=r'$A(\omega)$ from poles', color='red', linestyle='--')
        ax.plot(omega, sf.A_from_poles(omega).imag, label=r'$A(\omega)$ from poles', color='green', linestyle='--')
        ax.set_xlabel('Frequency (ω)')
        ax.set_ylabel(r'$A(\omega)$')


    # plot the poles of A(W)
    if(1):
        fig, ax = plt.subplots(1, 1, figsize=(6.5, 3.2))  # width, height in inches
        plt.rcParams.update({'font.size': 10})
        ax.scatter(sf.A_poles.real, sf.A_poles.imag, color='blue', marker='o', label='Poles of $A(\omega)$')
        ax.axhline(0, color='black', linewidth=0.5, linestyle='--')
        ax.axvline(0, color='black', linewidth=0.5, linestyle='--')
        ax.set_xlabel('Real Part')
        ax.set_ylabel('Imaginary Part')
        ax.set_title('Poles of $A(\\omega)$ in the Complex Plane')
        # label each pole with its residue
        for i in range(len(sf.A_poles)):
            ax.text(sf.A_poles[i].real, sf.A_poles[i].imag, f'{np.round(sf.A_res[i],2):.2e}', fontsize=8, ha='right')
        ax.legend()


    if(0):
        # plot S(w)
        fig, ax = plt.subplots(1, 1, figsize=(6.5, 3.2))  # width, height in inches
        omega = np.linspace(-300, 300, 10000)  # Frequency range
        ax.plot(omega, sf.S(omega).real, label=r'$S(\omega)$', color='blue', linestyle='-')
        ax.plot(omega, sf.S(omega).imag, label=r'$S(\omega)$', color='red', linestyle='--')
        ax.set_xlabel('Frequency (ω)')
        ax.set_ylabel(r'$S(\omega)$')


    plt.show()
    sys.exit()