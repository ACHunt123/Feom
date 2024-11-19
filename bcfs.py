from utils import orthotransmatgen
import numpy as np
import matplotlib.pyplot as plt


# J(w) parameters
eta = 0.1
gam = 0.1


# PLot the spectral density
if (0):
    w = np.linspace(0,2,1000)
    J = eta*gam*w/(w**2+gam**2)
    plt.plot(w,J)
    plt.show()

# Parameters for N-bead approximation and temperature
N=3
mu=N//2
beta=30
hbar=1

def NbeadTCF(N):
    # derived constants
    betaN = beta/N
    wN=1/(betaN*hbar)
    wk = np.array([2*wN*np.sin(np.pi*k/N) for k in range(0,mu+1)])

    # prefactors for the N-bead TCF
    d = np.zeros(mu+1)
    d0sum = np.sum(1/(gam**2*np.ones_like(wk) + wk**2))
    d[0] = eta/beta + 2*eta*gam**2/beta + d0sum
    d[1:] = -2*eta*gam/beta + wk[1:]/(gam**2*np.ones_like(wk[1:]) - wk[1:]**2)
    dI = -hbar*eta*gam/2

    # N-bead TCF
    t = np.linspace(0,10,1000)
    C = np.zeros_like(t,dtype=complex)

    C += (d[0] + dI*0+1.j)*np.exp(-gam*t)
    for k in range(1,mu+1):
        C += d[k]*np.exp(-wk[k]*t)

    return t,C


def matsTCF(N):
    # derived constants
    betaN = beta/N
    wN=1/(betaN*hbar)
    wn = np.array([2*wN*np.pi*k/N for k in range(0,mu+1)])

    # prefactors for the N-bead TCF
    d = np.zeros(mu+1)
    d0sum = np.sum(1/(gam**2*np.ones_like(wn) + wn**2))
    d[0] = eta/beta + 2*eta*gam**2/beta + d0sum
    d[1:] = -2*eta*gam/beta + wn[1:]/(gam**2*np.ones_like(wn[1:]) - wn[1:]**2)
    dI = -hbar*eta*gam/2

    # N-bead TCF
    t = np.linspace(0,10,1000)
    C = np.zeros_like(t,dtype=complex)

    C += (d[0] + dI*0+1.j)*np.exp(-gam*t)
    for k in range(1,mu+1):
        C += d[k]*np.exp(-wn[k]*t)

    return t,C


t,CN = NbeadTCF(N)
t,Cmats = matsTCF(N)

plt.plot(t,CN.real,label='N-bead')
plt.plot(t,Cmats.real,label='Matsubara')
plt.show()





