
##############################################
# Parser for my HEOM module
##############################################
import argparse
import numpy as np
from Feom.utils import out_filename,printparams
fs = 0.02418884254 #au
'''
===============================================================================================================
GENERAL PARAMETERS: 
---------------------------------------------------------------------------------------------------------------
hbar = 1
L           The depth of the ADO expansion
K           The number of exponential terms in the BCF
ns          Number of states to be propagated
beta        Inverse temperature
tmax        Maximum time to propagate
dt          Time step size
===============================================================================================================
POTENTIALS:
---------------------------------------------------------------------------------------------------------------
Spin boson:
Delta
eps
---------------------------------------------------------------------------------------------------------------
Harmonic oscillator:
omega
m
[dx]
[xmin]
[xmax]
===============================================================================================================
BATH PARAMETERS:
J(w) = (\pi/2) * \sum_{\alpha} \frac{c_\alpha^2}{m_\alpha \omega_\alpha} \delta(w - w_\alpha)
---------------------------------------------------------------------------------------------------------------
Debye bath:
J(w) = \frac{\eta\gamma\omega}{\omega^2 + \gamma^2}
eta         Coupling strength
gam         Cuttoff frequency
================================================================================================================
SWITCHES:
---------------------------------------------------------------------------------------------------------------
potname = ['harmonic','spinboson']              Potential name
bathname = ['debye','debyeCothpoles']           Bath name

bathmode = ['nbead','nmats','matsubara']        IF bathname == 'debye'
bathmode = ['Pade[N,N]','Pade[N-1,N]','AAA']    IF bathname == 'debyeCothpoles'

NOTE: nmats and nbead are if we chose n beads, then set the frequencies to be the Matsubara/RP frequencies
this means that the c0 coefifients will contain a finute sum. The matsubara option includes all poles,
giving the tan() in c0. This is the one used with Ishizki-Tanimura terminator.

LTCorr = ['NZ2','IT','PT2']    Whether to add Nakajima-Zwanwig, Ishizki-Tanimura or 2nd order terminator [not available for nbeads/nmats as there are no exra terms in BCF]

--print_ADOs        Print the ADOs to file every N (hardcoded) timesteps (if present, default False)
--noSIA             Use RK4 step instead of SIA step (used present, default False)
===============================================================================================================
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
    parser.add_argument("--LTCorr",type=str,default=argparse.SUPPRESS,help="Low temp correction")
    ### Switches
    parser.add_argument("--print_ADOs", action="store_true", help="Print the ADOs to file every N timesteps")
    parser.add_argument("--noSIA", action="store_true", help="Use RK4 step instead of SIA step")

    return parser.parse_args()

### Get all the parameters
params = parse_args()

### remove the unnecessary parameters
if params.potname == 'harmonic':
    del params.Delta
    del params.eps
elif params.potname == 'spinboson': 
    del params.omega
    del params.m
    del params.dx
    del params.xmin
    del params.xmax

### Override some of the parameters based on the bathname
if params.bathmode in ['nmats','nbead']: del params.LTCorr


### Derived parameters (not input)
params.nttot =int(params.tmax/params.dt)+1

