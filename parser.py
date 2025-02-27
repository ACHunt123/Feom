
##############################################
# Parser for my HEOM module
##############################################
import argparse,sys
import numpy as np
fs = 0.02418884254 #au


'''
Switches:
potname=['harmonic','spinboson'][1]
bathname = ['debye'][0]
bathmode = ['nbead','matsubara'][1]

Spin boson:
Delta
eps

Harmonic oscillator:
omega
m
[dx]
[xmin]
[xmax]

General: 
hbar=1
L = 3           # the depth of the ADO expansion
K = 3           #  the number of elements in the BCFs
ns = 2         # number of states to be propagated

Bath:
Delta=1
wc=1*Delta
eps=1*Delta
Lambda = 0.5*Delta
beta=0.25/Delta

eta = Lambda/2
gam = wc
'''
def parse_args():
    parser = argparse.ArgumentParser(description="Parser for Heom module")
    ### Potentials
    parser.add_argument("--potname",type=str,default='spinboson',help="Potential name")
    # Spin boson
    parser.add_argument("--Delta",type=float,default=1,help="Delta")
    parser.add_argument("--eps",type=float,default=1,help="eps")
    # Harmonic oscillator
    parser.add_argument("--omega",type=float,default=1,help="omega")
    parser.add_argument("--m",type=float,default=1741.1,help="m")
    parser.add_argument("--dx",type=float,default=0.01,help="dx")
    parser.add_argument("--xmin",type=float,default=-5,help="xmin")
    parser.add_argument("--xmax",type=float,default=5,help="xmax")
    ### General
    parser.add_argument("--hbar",type=float,default=1,help="hbar")
    parser.add_argument("--L",type=int,default=0,help="L")
    parser.add_argument("--K",type=int,default=0,help="K")
    parser.add_argument("--ns",type=int,default=10,help="number of states")
    parser.add_argument("--beta",type=float,default=0.25,help="beta")
    ### Bath
    parser.add_argument("--bathname",type=str,default='debye',help="Bath name")
    parser.add_argument("--bathmode",type=str,default='matsubara',help="Bath mode")
    parser.add_argument("--eta",type=float,default=0.25,help="bath strength")
    parser.add_argument("--gam",type=float,default=1,help="cuttoff frequency")
    return parser.parse_args()

params = parse_args()

### Propagation
params.tmax = 20
params.dt= 0.001
params.nttot =int(params.tmax/params.dt)+1
params.t_arr = np.arange(params.nttot)*params.dt

### more Parameters
    # beta600 = 526.2918822622093
    # beta300 = 1052.5837645244185
    # beta150 = 2105.167529048837

    # beta = beta150 #in Adam's code

    ### Bath parameters - Debye bath [from ADAM]
    # m =1741.1
    # d0 = 0.18748
    # alpha = 1.1605
    # diff = 2 * d0 *  alpha**2
    # omega = (diff/m)**(0.5)
    # eta_crit = 2*m*omega  #critical cutoff frequency 
    # eta_ADAM=2*eta_crit
    # eta = eta_ADAM*omega 
    # gam= omega

# Delta=1
# wc=1*Delta
# eps=1*Delta
# Lambda = 0.5*Delta
# beta=0.25/Delta

# eta = Lambda/2
# gam = wc

        # ### Hamiltonian and system setup - harmonic oscillator - same as Adam's code
        # m =1741.1
        # d0 = 0.18748
        # alpha = 1.1605
        # diff = 2 * d0 *  alpha**2
        # const = d0 *  alpha**2
        # omega = (diff/m)**(0.5)
