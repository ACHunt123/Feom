import numpy as np
from scipy.sparse import coo_matrix, eye, kron
from Feom.src.fermions.hierarchy.ADO_ops import generate_ado_raising_lowering_ops



def generate_liouvillian(sim):

    ### Precalculate all of the matrices etc
    # parameters
    # length for the normally used terms
    K = len(sim.bath.C_ks_plus)
    # extra length if RWA added
    if hasattr(sim.bath, "gam_js"):
        J = len(sim.bath.gam_js)
    else:
        J=0

    L = sim.params.L
    # Raising/Lowering and number operators in the ADO space 
    ADO_ops,indices_dict = generate_ado_raising_lowering_ops(2*K+J,L)
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

    # 'Normal' treatment of C_ks modes (two channels for each)
    for k in range(K):
        # The ordering of the + - modes in ADO space
        kap_plus=k*2      #kappa plus
        kap_mnus=k*2+1    #kappa minus

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
        L_terms.append(-1.j*kron(P_modes[kap_plus]@A[kap_plus], VmL, format='coo'))
        L_terms.append(1.j*kron(P_modes[kap_plus]@P_global@A[kap_plus], VmR, format='coo'))
        # s=-
        L_terms.append(-1.j*kron(P_modes[kap_mnus]@A[kap_mnus], VpL, format='coo'))
        L_terms.append(1.j*kron(P_modes[kap_mnus]@P_global@A[kap_mnus], VpR, format='coo'))


        ## Add lowering terms to the Liouvillian
        # s=+
        L_terms.append(kron(-1.j*Cp*P_modes[kap_plus]@Adag[kap_plus], VpL, format='coo'))
        L_terms.append(kron(-1.j*Cm.conj()*P_modes[kap_plus]@P_global@Adag[kap_plus], VpR, format='coo'))
        # s=-
        L_terms.append(kron(-1.j*Cm*P_modes[kap_mnus]@Adag[kap_mnus], VmL, format='coo'))
        L_terms.append(kron(-1.j*Cp.conj()*P_modes[kap_mnus]@P_global@Adag[kap_mnus], VmR, format='coo'))
        
        
        ## Add on the diagonal (-n gamma) damping
        L_terms.append(-gp*kron(AdagA[kap_plus], I_sys, format='coo'))  
        L_terms.append(-gm*kron(AdagA[kap_mnus], I_sys, format='coo'))

    # RWA treatment of C_js terms (one channel for each)
    # see the non-markovian quantum skin effect paper for RWA example paper for bosons.
    if not sim.params.IT_RWAterms: # Do the anti-RWA approximation on the j terms
        for j in range(J):
            kap = 2*K+j

            # get the coupling operators (p,m=+,- and L/R for left/right acting)
            Vp= sim.pot.V_js_plus[j] 
            VpL= kron(Vp,I, format='coo')
            VpR= kron(I,Vp.T, format='coo')
            Vm= sim.pot.V_js_mnus[j] 
            VmL= kron(Vm,I, format='coo')
            VmR= kron(I,Vm.T, format='coo')

            # get the coefficients (abbreviated)
            Cp=np.complex128(sim.bath.C_js_plus[j])
            Cm=np.complex128(sim.bath.C_js_mnus[j])
            g=np.complex128(sim.bath.gam_js[j])


            ### Add on the off-diagonal (nk+ and nk-) coupling
            ## Add raising terms to the Liouvillian
            L_terms.append(-1.j*kron(P_modes[kap]@A[kap],(VmL+VpL), format='coo'))
            L_terms.append(1.j*kron(P_modes[kap]@P_global@A[kap], (VmR+VpR), format='coo'))

            ## Add lowering terms to the Liouvillian
            L_terms.append(kron(-1.j*P_modes[kap]@Adag[kap], (Cm*VmL+Cp*VpL), format='coo'))
            L_terms.append(kron(-1.j*P_modes[kap]@P_global@Adag[kap], (Cp.conj()*VmL+Cm.conj()*VpL), format='coo'))
            
            ## Add on the diagonal (-n gamma) damping
            L_terms.append(-g*kron(AdagA[kap], I_sys, format='coo'))  

        if J!=0 and np.all(np.array(sim.bath.C_js_plus)!=-np.array(sim.bath.C_js_mnus).conj()):
            print('Need to add on the IT closure, as terms do not cancel')
            exit()

    # Markovian truncation for C_js (no channel for each)    
    else: # Do ishizaki-Tanimura on the rotating wave terms instead
        for j in range(J):

            # get the coupling operators (p,m=+,- and L/R for left/right acting)
            Vp= sim.pot.V_js_plus[j] 
            VpL= kron(Vp,I, format='coo')
            VpR= kron(I,Vp.T, format='coo')
            Vm= sim.pot.V_js_mnus[j] 
            VmL= kron(Vm,I, format='coo')
            VmR= kron(I,Vm.T, format='coo')

            # get the coefficients (abbreviated)
            Cp=sim.bath.C_js_plus[j]
            Cm=sim.bath.C_js_mnus[j]
            g=sim.bath.gam_js[j]
            
            ### Add on the off-diagonal (nk+ and nk-) coupling
            ## Add raising terms to the Liouvillian
            # s=+
            L_terms.append(-(1./g)*kron(I_ado, (Cp*VmL@VpL-Cm.conj()*VpR@VmR), format='coo'))
            L_terms.append((1./g)*kron(P_global,(Cm.conj()*VmL@VpR + Cp*VpL@VmR), format='coo'))
            # s=-
            L_terms.append(-(1./g)*kron(I_ado, (Cm*VpL@VmL-Cp.conj()*VmR@VpR), format='coo'))
            L_terms.append((1./g)*kron(P_global,(Cp.conj()*VpL@VmR + Cm*VmL@VpR), format='coo'))

    # set the outputs
    L_combined = sum(L_terms)
    # Now convert to the correct size CSR matrix for the solver
    sim.Liouvillian = L_combined.tocsr()
    sim.params.Nados=Nados
    sim.params.Ntot=sim.params.Nados*sim.params.ns**2
    
