import numpy as np
from scipy.sparse import coo_matrix, eye, kron
from Feom.src.hierarchy.ADO_ops import generate_ado_raising_lowering_ops

def organize_exponents(gam_ks,tol=1e-6):
    """
    This function categorizes BCF exponents (gam_ks) based on their imaginary parts. 
    Real exponents (e.g., Matsubara terms) are identified as purely exponential 
    decay. Complex exponents are paired by matching equal real parts and 
    opposite imaginary parts to facilitate the construction of conjugate 
    coupling operators.

    Args:
        gam_ks (array_like): Complex exponents (decay rates) of the bath 
            correlation function.
        tol (float, optional): Absolute tolerance for determining if the 
            imaginary part is zero or if real parts match. Defaults to 1e-6.

    Returns:
        k_decay_inds (numpy.ndarray): Indices of purely real exponents 
            representing simple exponential decay.
        k_oscil_inds (numpy.ndarray): Nx2 array where each row contains the 
            index pair [i, j] such that gam_ks[i] and gam_ks[j] are 
            complex conjugates.

    Raises:
        ValueError: If the exponents cannot be uniquely partitioned into 
            real terms and conjugate pairs, or if indices are duplicated.
    """
    ## get the masks for the pure exponential decay and the +- damped exponentials
    k_decay_inds = np.flatnonzero(np.abs(np.imag(gam_ks)) <= tol)
    k_plus_inds  = np.flatnonzero(np.imag(gam_ks) >= tol)
    k_minus_inds = np.flatnonzero(np.imag(gam_ks) <= -tol)

    # exit early if all of the exponentials are pure damped
    if np.array_equal(np.sort(np.concatenate([k_decay_inds])),np.arange(len(gam_ks))):
        return k_decay_inds, np.empty((0, 2), dtype=int)
    # check that the masks are containing all of the elements
    if not np.array_equal(np.sort(np.concatenate([k_decay_inds,k_plus_inds,k_minus_inds])),np.arange(len(gam_ks))):
        raise ValueError("Index masks do not cover all poles exactly once.")
    
    ## match the real parts of the +- damped exponentials
    k_oscil_inds=[]
    for k_plus_ind_i in k_plus_inds:
        for k_minus_ind_i in k_minus_inds:
            # get the gam_ks corresponding to these indices
            gam_ki_plus = gam_ks[k_plus_ind_i]
            gam_ki_minus = gam_ks[k_minus_ind_i]
            # check if the real parts are the same (then they are pairs)
            if np.real(gam_ki_plus)==np.real(gam_ki_minus):
                # add the pair of indices to the list
                k_oscil_inds.append([k_plus_ind_i,k_minus_ind_i])
    k_oscil_inds=np.asarray(k_oscil_inds,dtype=int)
    # Check they all got paired up
    if not np.array_equal(np.sort(np.concatenate([k_decay_inds,k_oscil_inds[:,0],k_oscil_inds[:,1]])),np.arange(len(gam_ks))):
        raise ValueError("Index masks do not cover all poles exactly once.")

    return k_decay_inds, k_oscil_inds



def generate_liouvillian(sim):

    ### Precalculate all of the matrices etc
    # parameters
    K = len(sim.bath.C_ks)
    L = sim.params.L
    # Raising/Lowering and number operators in the ADO space
    ADO_ops,indices_dict = generate_ado_raising_lowering_ops(K,L)
    A = ADO_ops.A
    Adag = ADO_ops.Adag
    AdagA = ADO_ops.AdagA
    Nados = len(indices_dict)
    # Identity matrices
    I_ado = eye(Nados, format='coo', dtype=np.complex128)
    I_sys = eye(sim.params.ns**2, format='coo', dtype=np.complex128) # LIOUVILLE 
    # System Liouville space operators
    I= np.eye(sim.pot.ns) #the hilbert space identity (NOT I_sys)
    H= sim.pot.H_mat #NOTE THIS IS THE RENORMALIZED HAMILTONIAN

    ### Build the Liouvillian

    ## Diagonal, ado-index-independent terms
    L_terms = [] # This list will hold all components of the Liouvillian

    # add system free liouviliian matrix
    Lsys0 = coo_matrix(-1.j*(np.kron(H,I) - np.kron(I,H.T)), dtype=np.complex128)
    L_terms.append(kron(I_ado, Lsys0, format='coo'))
    # add on terminator
    Xi= sim.Xi
    if sim.params.LTCorr == 'same_for_each_ADO':
        L_terms.append(kron(I_ado, Xi, format='coo'))
    elif sim.params.LTCorr == 'different_for_each_ADO':
        print('not implemented yet')

    ## ado-index-dependent terms
    # get the indices of the pure exponential/oscillitory pairs of bath exponentials
    k_decay_inds, k_oscil_inds = organize_exponents(sim.bath.gam_ks)
    ## pure decay modes first
    for ki in k_decay_inds:
        # get the coupling operators
        V= sim.pot.V_ks[ki] 
        VL= kron(V,I, format='coo')
        VR= kron(I,V.T, format='coo')
        Vx= VL-VR

        # get the coefficients
        Ck=sim.bath.C_ks[ki]
        gamk=sim.bath.gam_ks[ki]
        # Add on the off-diagonal (nk+ and nk-) coupling
        Lk_dn = -(1.j/np.sqrt(np.abs(Ck)))*(Ck*VL-Ck.conj()*VR) #coupling to ados lower in excitation
        Lk_up = -1.j*np.sqrt(np.abs(Ck))*Vx                     #coupling to ados higher in excitation
        L_terms.append(kron(Adag[ki], Lk_dn, format='coo'))
        L_terms.append(kron(A[ki], Lk_up, format='coo'))
        # Add on the diagonal (-n gamma) damping
        L_terms.append(-gamk*kron(AdagA[ki], I_sys, format='coo'))  

    ## oscilitory modes now
    # we assume sum-difference expansion from of the IF expansion
    # generalization to the forward/backward path separation is
    # a simple rotation in creation annihilation operators of ADOs
    # and is a simple couple extra lines of code
    for k_plus,k_minus in k_oscil_inds:
        # get the coupling operators + check if they are the same for k_plus and k_minus
        V= sim.pot.V_ks[k_plus] 
        if not np.array_equal(V,sim.pot.V_ks[k_minus]):
            raise ValueError(f"Coupling operators at indices {k_plus} and {k_minus} must be identical.")
        VL= kron(V,I, format='coo')
        VR= kron(I,V.T, format='coo')
        Vx= VL-VR
    
        # get the coefficients
        Ck_plus=sim.bath.C_ks[k_plus]
        Ck_minus=sim.bath.C_ks[k_minus]
        gamk_plus=sim.bath.gam_ks[k_plus]
        gamk_minus=sim.bath.gam_ks[k_minus]
        
        ## Add on the off-diagonal (nk+ and nk-) coupling
        # Calculate the average weighting for the pairs of modes
        Wk=(np.sqrt(np.abs(Ck_plus))+np.sqrt(np.abs(Ck_minus)))/2.
        # Lk_up_plus = -1.j * np.sqrt(np.abs(Ck_plus)) * Vx
        # Lk_up_minus = -1.j * np.sqrt(np.abs(Ck_minus)) * Vx
        # Lk_dn_plus = -(1.j / np.sqrt(np.abs(Ck_plus))) * (Ck_plus * VL - np.conj(Ck_minus) * VR)
        # Lk_dn_minus = -(1.j / np.sqrt(np.abs(Ck_minus))) * (Ck_minus * VL - np.conj(Ck_plus) * VR)
        Lk_up_plus = -1.j * Wk * Vx
        Lk_up_minus = -1.j * Wk * Vx
        Lk_dn_plus = -(1.j / Wk) * (Ck_plus * VL - np.conj(Ck_minus) * VR)
        Lk_dn_minus = -(1.j / Wk) * (Ck_minus * VL - np.conj(Ck_plus) * VR)
        # Add raising terms to the Liouvillian
        L_terms.append(kron(A[k_plus], Lk_up_plus, format='coo'))
        L_terms.append(kron(A[k_minus], Lk_up_minus, format='coo'))
        # Add lowering terms to the Liouvillian
        L_terms.append(kron(Adag[k_plus], Lk_dn_plus, format='coo'))
        L_terms.append(kron(Adag[k_minus], Lk_dn_minus, format='coo'))

        ## Add on the diagonal (-n gamma) damping
        L_terms.append(-gamk_plus*kron(AdagA[k_plus], I_sys, format='coo'))  
        L_terms.append(-gamk_minus*kron(AdagA[k_minus], I_sys, format='coo'))

    # set the outputs
    L_combined = sum(L_terms)
    # Now convert to the correct size CSR matrix for the solver
    sim.Liouvillian = L_combined.tocsr()
    sim.params.Nados=Nados
    sim.params.Ntot=sim.params.Nados*sim.params.ns**2
    
