import numpy as np
from scipy.sparse import coo_matrix, eye, kron
from Feom.src.fermions.hierarchy.ADO_ops import generate_ado_raising_lowering_ops



def generate_liouvillian(sim):

    ### Precalculate all of the matrices etc
    # parameters
    K = len(sim.bath.C_ks_plus)+len(sim.bath.C_ks_mnus)
    L = sim.params.L
    # Raising/Lowering and number operators in the ADO space 
    ADO_ops,indices_dict = generate_ado_raising_lowering_ops(K,L)
    A = ADO_ops.A
    Adag = ADO_ops.Adag
    AdagA = ADO_ops.AdagA
    P_global=ADO_ops.P_global   # Global permutational sign (-1)^n where n is the tier of the ADO
    P_modes=ADO_ops.P_modes     # Local sign (-1)^(n-k) where n is the tier, and k is for the kth excitation
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
    # Xi= sim.Xi
    # if sim.params.LTCorr == 'same_for_each_ADO':
    #     L_terms.append(kron(I_ado, Xi, format='coo'))
    # elif sim.params.LTCorr == 'different_for_each_ADO':
    #     print('not implemented yet')

 
    for k in range(K//2):
        # The ordering of the + - modes in ADO space
        k_plus=k*2
        k_mnus=k*2+1

        # get the coupling operators (p,m=+,- and L/R for left/right acting)
        Vp= sim.pot.V_ks_plus[k] 
        VpL= kron(Vp,I, format='coo')
        VpR= kron(I,Vp.T, format='coo')
        Vm= sim.pot.V_ks_mnus[k] 
        VmL= kron(Vm,I, format='coo')
        VmR= kron(I,Vm.T, format='coo')

        # get the coefficients (abbreviated)
        Cp=sim.bath.C_ks_plus[k]
        Cm=sim.bath.C_ks_mnus[k]
        gp=sim.bath.gam_ks_plus[k]
        gm=sim.bath.gam_ks_mnus[k]

        
        ### Add on the off-diagonal (nk+ and nk-) coupling
        ## Add raising terms to the Liouvillian
        # s=+
        L_terms.append(-1.j*kron(P_modes[k_plus]@A[k_plus], VmL, format='coo'))
        L_terms.append(1.j*kron(P_modes[k_plus]@P_global@A[k_plus], VmR, format='coo'))
        # s=-
        L_terms.append(-1.j*kron(P_modes[k_mnus]@A[k_mnus], VpL, format='coo'))
        L_terms.append(1.j*kron(P_modes[k_mnus]@P_global@A[k_mnus], VpR, format='coo'))


        ## Add lowering terms to the Liouvillian
        # s=+
        L_terms.append(kron(-1.j*P_modes[k_plus]@Adag[k_plus], Cp*VpL, format='coo'))
        L_terms.append(kron(-1.j*P_modes[k_plus]@P_global@Adag[k_plus], Cm.conj()*VpR, format='coo'))
        # s=-
        L_terms.append(kron(-1.j*P_modes[k_mnus]@Adag[k_mnus], Cm*VmL, format='coo'))
        L_terms.append(kron(-1.j*P_modes[k_mnus]@P_global@Adag[k_mnus], Cp.conj()*VmR, format='coo'))
        
        
        ## Add on the diagonal (-n gamma) damping
        L_terms.append(-gp*kron(AdagA[k_plus], I_sys, format='coo'))  
        L_terms.append(-gm*kron(AdagA[k_mnus], I_sys, format='coo'))

    # set the outputs
    L_combined = sum(L_terms)
    # Now convert to the correct size CSR matrix for the solver
    sim.Liouvillian = L_combined.tocsr()
    sim.params.Nados=Nados
    sim.params.Ntot=sim.params.Nados*sim.params.ns**2
    
