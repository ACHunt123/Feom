#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from Feom.src.initialize.setup import Setup,SimConfig
from A34BCFs import A4_BCF,A3_BCF,make_bcf_func, Sbeta_exact

# plotting switches
plot_fits=1
# general params
L=4
W=6
ns=2
beta=1000
tmax= 10
dt=0.01
# system params
Delta=2
eps=1
# bath params
eta_DL=2
gam_DL=1
# fit parameters
w_max=500
w_max=200
N_support=100000


### Setup the system
H_mat = np.zeros((ns,ns),dtype=complex)
H_mat[0,0] = -eps
H_mat[1,1] = eps
H_mat[1,0] = Delta
H_mat[0,1] = Delta
# perturbation matrix
q_mat = np.zeros((ns,ns),dtype=complex)
q_mat[0,0] = -1
q_mat[1,1] = 1

### Setup the bath
if(1): # overdamped dl
    Jw_pos_residues = [eta_DL*gam_DL/2]
    Jw_pos_poles=[1.j*gam_DL]
else: # try underdamped dl
    Omega_D = 2
    gamma_D = 1
    lambda_D = 1
    Xi_D = np.sqrt(Omega_D**2 - gamma_D**2)
    Jw_pos_poles = [Xi_D + 1.j * gamma_D,  -Xi_D + 1.j * gamma_D   ]
    res_val = (0.5* lambda_D * Omega_D**2) / (1.j * Xi_D)
    Jw_pos_residues= [res_val, -res_val]


if(0): #plot the different approximations of the correlation functions
    C_ks,gam_ks,zeta= A4_BCF(Jw_pos_poles, Jw_pos_residues,W,beta,doplot=plot_fits)
    print(f'mength {len(C_ks)}')
    A4_approx=make_bcf_func(C_ks,gam_ks)

    C_ks,gam_ks,zeta= A3_BCF(Jw_pos_poles, Jw_pos_residues,W,beta,doplot=plot_fits)
    print(f'mength {len(C_ks)}')

    A3_approx=make_bcf_func(C_ks,gam_ks)

    import matplotlib.pyplot as plt
    t_arr=np.linspace(0,100,1000)
    # plot the approximations
    plt.plot(t_arr,A4_approx(t_arr).real,label='A4')
    plt.plot(t_arr,A3_approx(t_arr).real,label='A3')
    # plot the exact (needs a FT)
    Sbeta,omega=Sbeta_exact(Jw_pos_poles, Jw_pos_residues,beta,w_max=w_max,N_support=N_support)
    dw=omega[2]-omega[1]
    C_exact=np.sum(Sbeta[:,None]*np.exp(1j*omega[:,None]*t_arr[None,:]),axis=0)*dw
    plt.plot(t_arr,C_exact.real,label='exact',ls='--',zorder=100)

    plt.legend()
    plt.show()
    exit()

# calc the params for both
A3params = (A3_BCF(Jw_pos_poles, Jw_pos_residues,W,beta,doplot=plot_fits,w_max=w_max,N_support=N_support))
A4params = (A4_BCF(Jw_pos_poles, Jw_pos_residues,W,beta,doplot=plot_fits,w_max=w_max,N_support=N_support))

#run both
plt.cla()
for (C_ks,gam_ks,zeta),name in zip([A3params,A4params],['A3','A4']):
    
    ### Setup the terminator
    I = np.eye(ns)
    Vcross = np.kron(q_mat,I) - np.kron(I,q_mat.T)  # commutator superoperator for the system-bath coupling operator
    Xi= -1 * (zeta/2) * Vcross @ Vcross # Add on the terminator contribution from the delta function in BCF
    ### Build the dictionaries
    sys_args = {'s_mat': q_mat,'H_mat': H_mat,}
    bath_args = {'C_ks':C_ks,'gam_ks':gam_ks,'zeta':zeta}
    params_args = {'dt': dt,'tmax': tmax,'L': L}
    terminator_args = {'correction_type': 'same_for_each_ADO','Xi':Xi}
    ### Initialize
    my_config=SimConfig(sys_args, bath_args, params_args,terminator_args)
    sim = Setup(my_config)
    # set initial conditions
    rhos_0 = np.zeros((2, 2), dtype=complex)
    rhos_0[1, 1] = 1.0 # Population in state 1
    sim.set_initial_ADOs(rhos_0,'0th')
    ### Run
    sim.generate_input_files()
    sim.insert_executable()
    sim.go(cleanup=1)
    # Process results
    rho11 = sim.rho[:,1,1]  
    rho00 = sim.rho[:,0,0]
    rho10 = sim.rho[:,1,0]
    rho01 = sim.rho[:,0,1]

    processed_data = np.zeros((len(sim.t_arr),5),dtype=complex)
    processed_data[:,0] = sim.t_arr
    processed_data[:,1] = rho11
    # processed_data[:,1] = rho11 - rho00  # <s_z>
    # processed_data[:,2] = 1.j*(rho10 - rho01)  # <s_y>
    # processed_data[:,3] = (rho10 + rho01)  # <s_x>   
    # processed_data[:,4] = rho11 # Site 1 population
    # data_labels = '\n,Time /a.u. <s_z> <s_y> <s_x> P11'
    data_labels = '\n,Time /a.u.  rho11 - - - '
    # np.savetxt('exact(A4_W=22).dat',processed_data.real,header=data_labels)

    plt.plot(sim.t_arr,rho11,label=f'{name} W={len(sim.bath.C_ks)}')

exact=np.loadtxt('exact(A4_W=22).dat')
plt.title(f'Comparison between AAA and A4 with L={L}')
plt.plot(exact[:,0],exact[:,1],label='exact(A4 with W=20)',zorder=100, ls='--')
plt.legend()

plt.ylim(0,1)
plt.xlim(0,tmax)
plt.savefig('plot.pdf')
plt.show()



