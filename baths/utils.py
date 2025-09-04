import numpy as np
np.set_printoptions(precision=5,linewidth=200,suppress=True)
from types import SimpleNamespace
import Feom.baths as baths
import sys,copy

def get_C_UDs(params):
    ''' Generate the coefficients C_U, c_D_LEFT, c_D_RIGHT from the C_ks for the bath (that are used in the FEOM code)'''
    params.c_U = np.zeros((params.N_exp,params.L+1),dtype=complex)
    params.c_D_LEFT = np.zeros((params.N_exp,params.L+1),dtype=complex)
    params.c_D_RIGHT = np.zeros((params.N_exp,params.L+1),dtype=complex)
    for ki in range(params.N_exp):
        for nk in range(params.L+1):
            params.c_U[ki,nk] = np.sqrt((nk+1)*abs(params.C_ks[ki]))
            if abs(params.C_ks[ki]) < 1e-10:
                params.c_D_LEFT[ki,nk] = 0.0
                params.c_D_RIGHT[ki,nk] = 0.0
                print(f'Warning: C_ks({ki}) is zero, setting superoperator terms to zero')
            else:
                params.c_D_LEFT[ki,nk] = -np.sqrt(nk/abs(params.C_ks[ki]))*params.C_ks[ki]
                params.c_D_RIGHT[ki,nk] = np.sqrt(nk/abs(params.C_ks[ki]))*np.conj(params.C_ks[ki])
    return 


def generate_Terminator(sim):
    '''Function to generate terminator for HEOM.
    IT - ishizaki-Tanimura terminator (same for each ADO)
    PT2 - 2nd order perturbative terminator (same for each ADO)
    NZ2 - Nakajima-Zwanwig terminator (different for each ADO, more expensive)


    self.init_lowtcoef is the low-temp correction either due to IT, or the k term from AAA and Pade[N/N] it will be added to
    It is a class, as it need to store its type and size

    C_exact(t) = C_approx(t) + C_mats_inf(t) - C_approx(t),
    where C_approx(t) is the approximation we are using (e.g. Pade etc)
    then we expand
    C_exact(t) = C_approx(t) + [C_mats_Kbig(t) - C_approx(t)] + C_mats_K>Kbig(t)
                                <------- deltaC(t) ------->             
    
    deltaC(t)           : treated peturbatively (PT2, NJZ) or markovianly (IT)
    C_mats_K>Kbig(t)    : treated Markovianly (the IT terminator)
    C_approx(t)         : explicitly propagated in the HEOM
    '''
    params = sim.params
    pot= sim.pot
    bath = sim.bath

    ### Get termination coefficients (the FAY Way)
    def get_termination_coefs(sim):
        ''' Get the Cks and gamks for the terminated frequencies corresponding to 
        deltaC(t) = C_{mats Kbig}(t) - C_{approx}(t)
        '''
        # setup the matsubara bath object C_{mats Kbig}(t)
        matsbath_params=copy.deepcopy(sim.params)
        matsbath_params.bathmode='matsubara'
        Kbig=500        # up to K=500, 
        matsbath_params.K=Kbig
        matsbath = baths.getbath('debye')(matsbath_params)
        
        # Get the Cks and gamks for the deltaC(t)
        if bath.mode != 'matsubara':  
            delta_Cks = np.append(-sim.bath.C_ks[:sim.bath.N_exp_prop],matsbath.C_ks)
            delta_gamks = np.append(sim.bath.gam_ks[:sim.bath.N_exp_prop],matsbath.gam_ks)
        else:
            delta_Cks = matsbath.C_ks[sim.bath.N_exp_prop:]
            delta_gamks = matsbath.gam_ks[sim.bath.N_exp_prop:]

        # Get the Markovian terms from C_mats_K>Kbig(t) for the IT terminator
        mark_corr = matsbath.lowTcoef
        print(f'Low temp correction term from terminated frequencies: {mark_corr}')

        return delta_Cks, delta_gamks, mark_corr



    ### Get the Liouvillian and transformation matrices (for the NZ2 terminator)
    def getL(params,pot):
        ''' Get the Liouvillian matrix for the system Hamiltonian, and transformation matrices'''
        I = np.eye(params.ns)
        Lsys = -1.j*(np.kron(pot.H_mat,I) - np.kron(I,pot.H_mat.T))/params.hbar
        eigvals, eigvecs = np.linalg.eig(Lsys)
        Lams = np.diag(eigvals)
        Pis = eigvecs
        Pis_inv = np.linalg.inv(Pis)
        if(0):#test the eigen-decomposition
            L_reconstructed = Pis @ Lams @ Pis_inv
            error = np.linalg.norm(Lsys - L_reconstructed)
            print(f"Eigen-decomposition reconstruction error: {error:.2e}")
        return SimpleNamespace(Lsys=Lsys,Lams=np.diag(eigvals),Pis=eigvecs,Pis_inv=np.linalg.inv(eigvecs))



    def get_Xi_n(sim,I,Ldata):
        ''' Get the Xi_n matrix for a given ADO I (list of indices)
        for the terminated frequencies'''

        # first get the sum of the gam_ks*n_k for this ADO
        nks = sim.ADO_index[I,:]
        gamks = bath.gam_ks[:bath.N_exp_prop]
        gamma_n = np.sum(nks*gamks)  # the sum of n_k * gam_k
        # now calculate the terminator
        I_s_hilbert = np.eye(params.ns)
        VL = np.kron(pot.s_mat,I_s_hilbert)
        VR = np.kron(I_s_hilbert,pot.s_mat.T)
        V_x = VL-VR

        # calculate the sum
        Xi_n=np.zeros_like(Ldata.Lsys)
        for Ck,gamk in zip(bath.Cks_term,bath.gamks_term):
            if abs(Ck) < 1e-10:
                print('Warning: Ck is zero, skipping term in terminator')
                continue
            L_k_plus = -1.j*(Ck*VL -Ck.conj()*VR)/(params.hbar)
            L_k_minus= -1.j*V_x/params.hbar
            fraction = np.diag(1/(gamk + gamma_n - np.diag(Ldata.Lams)))
            Xi_n += L_k_minus @ Ldata.Pis @ (fraction) @ Ldata.Pis_inv @ L_k_plus

        # remove the k term contribution if it exists (it is added separately)
        # this is because dC(t) will have the k term in it too, so we need to remove it from what Xi would be if no terminator was added
        k = getattr(sim.bath, 'k', 0)                                           # get the constant term in the BCF if it exists, otherwise 0
        k_term_coef = -2*sim.bath.eta*sim.bath.gam*k/(params.beta*params.hbar**2)   # the constant term in the BCF gives a low temperature correction term
        Xi_n -= k_term_coef * V_x @ V_x
        return Xi_n

    def get_IT(params,pot,bath):
        ''' Get the Ishizaki-Tanimura terminator for the terminated frequencies
        '''
        I = np.eye(params.ns)
        Vcross = np.kron(pot.s_mat,I) - np.kron(I,pot.s_mat.T)  # commutator superoperator for the system-bath coupling operator
        IT_term_coef =  -np.sum(bath.Cks_term/bath.gamks_term)/params.hbar**2   # the IT terminator coefficient
        return IT_term_coef * Vcross @ Vcross


    ### Calculate the terminator
    Xi = np.zeros((1,params.ns**2,params.ns**2),dtype=complex)
    I = np.eye(params.ns)
    Vcross = np.kron(pot.s_mat,I) - np.kron(I,pot.s_mat.T)  # commutator superoperator for the system-bath coupling operator
    # Add on the terminator contribution from constant term k in BCF
    k = getattr(bath, 'k', 0)                                           # get the constant term in the BCF if it exists, otherwise 0
    k_term_coef = -2*bath.eta*bath.gam*k/(params.beta*params.hbar**2)   # the constant term in the BCF gives a low temperature correction term
    Xi += k_term_coef * Vcross @ Vcross

    ### Add on the terminator contribution from the terminated frequencies
    bath.LTCorr = getattr(params, 'LTCorr', None)
    # get the Cks and gam_ks for which we have terminated
    bath.Cks_term,bath.gamks_term,bath.mark_corr = get_termination_coefs(sim)

    if bath.LTCorr is None:
        pass   

    elif bath.LTCorr == 'PT2':   # Second order perturbative terminator on the terminated frequencies
        Ldata = getL(params,pot)
        Xi += get_Xi_n(sim,0,Ldata)                 # add Xi0 of tom's for each ADO
        Xi += bath.mark_corr * Vcross @ Vcross      # add markovian term for the frequencies higher than those included in the PT2 terminator

    elif bath.LTCorr == 'IT':    # Markovian (Ishizaki-Tanimura) terminator on the terminated frequencies
        Xi += get_IT(params,pot,bath)               # Calculate the IT terminator from the terminated frequencies
        Xi += bath.mark_corr * Vcross @ Vcross      # add markovian term for the frequencies higher than those included in the IT terminator

    elif bath.LTCorr == 'NZ2':   # Tom Fay's Nakajima-Zwanwig terminator on the terminated frequencies
        Xi_n = np.zeros((params.Imax,params.ns**2,params.ns**2),dtype=complex)  # one xi for EACH ADO (will overwrite Xi)
        Ldata = getL(params,pot)
        for I in range(params.Imax):
            Xi_n[I,:,:] += get_Xi_n(sim,I,Ldata)
            Xi_n[I,:,:] += Xi[0,:,:]                         # add on the contributions from the constant term k onto each Xi_n
            Xi_n[I,:,:] += bath.mark_corr * Vcross @ Vcross  # add markovian term for the frequencies higher than those included in the NZ2 terminator
            #NOTE: FAY does not include the Markovian term in his NZ2 terminator, but it should be there. see benchmark results (remove addition for exact agreement)
        # override Xi to be the full set of ADO-specific terminators
        Xi = Xi_n   

    else:
        raise ValueError('Invalid LTCorr type')
    
    ### Assign a dummy variable if the terminator is zero
    if (abs(Xi)==0).all():
        sim.Xi = np.zeros((1,1,1),dtype=complex) # this is just a dummy var
    else:
        sim.Xi = Xi
    return






