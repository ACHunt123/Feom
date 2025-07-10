
##############################################
# Parser for my HEOM module
##############################################
import argparse,sys
import numpy as np
from Feom.utils import out_filename,printparams
fs = 0.02418884254 #au


'''
Switches:
potname=['harmonic','spinboson'][1]
bathname = ['debye','debye'][0]

bathmode = ['nbead','matsubara'] IF bathname == 'debye'
bathmode = ['Pade[N,N]','Pade[N-1,N]','AAA'] IF bathname == 'debyeCothpoles'

lowTCorr = 0 # 1 to add low temperature corrections, 0 to not add them [overridden for pade and AAA modes]

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
L = 3           The depth of the ADO expansion
K = 3           The number of Matsubara terms in the BCF/ Total number of modes for AAA expansion (must be even)
ns = 2          Number of states to be propagated

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
    ### General
    parser.add_argument("--L",type=int,default=0,help="L")
    parser.add_argument("--K",type=int,default=0,help="K")
    parser.add_argument("--ns",type=int,default=10,help="number of states")
    parser.add_argument("--beta",type=float,default=0.25,help="beta")
    parser.add_argument("--hbar",type=float,default=1,help="hbar")
    parser.add_argument("--tmax",type=float,default=20,help="tmax")
    parser.add_argument("--dt",type=float,default=0.001,help="dt")
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
    ### Bath
    parser.add_argument("--bathname",type=str,default='debye',help="Bath name")
    parser.add_argument("--bathmode",type=str,default='matsubara',help="Bath mode")
    parser.add_argument("--eta",type=float,default=0.25,help="bath strength")
    parser.add_argument("--gam",type=float,default=1,help="cuttoff frequency")
    ### Switches
    parser.add_argument("--lowTCorr",type=int,default=0,help="add low temp corrections?") # NEED TO CHANGE THIS TO BOOLEAN
    # parser.add_argument("--lowTCorr", action="store_true", help="Add low temperature corrections if this flag is set.")
    parser.add_argument("--print_ADOs", action="store_true", help="Print the ADOs to file every N timesteps")
    parser.add_argument("--prune", action="store_true", help="Prune the ADOs dynamically during propagation")

    return parser.parse_args()

### Get all the parameters
params = parse_args()
### remove the unnecessary ones
if params.potname == 'harmonic':
    del params.Delta
    del params.eps
elif params.potname == 'spinboson': 
    del params.omega
    del params.m
    del params.dx
    del params.xmin
    del params.xmax

### Write the parameters to a file and filename
params.header = printparams(params)
params.out_name = out_filename(params)


### Derived parameters (not input)
params.nttot =int(params.tmax/params.dt)+1
params.t_arr = np.arange(params.nttot)*params.dt

