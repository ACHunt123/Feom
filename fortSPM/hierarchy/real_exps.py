import numpy as np
from scipy.sparse import csr_matrix, eye, kron
from Feom.fortSPM.hierarchy.ADO_ops import generate_ado_raising_lowering_ops


def generate_liouvillian(sim):

    K = len(sim.bath.C_ks)
    L = sim.params.L

    # Raising/Lowering and number operators in the ADO space
    ADO_ops,indices_dict = generate_ado_raising_lowering_ops(K,L)
    A = ADO_ops.A
    Adag = ADO_ops.Adag
    AdagA = ADO_ops.AdagA
    Nados = len(indices_dict)

    # Identity matrices
    I_ado = eye(Nados, format='csr', dtype=np.complex128)
    I_sys = eye(sim.params.ns**2, format='csr', dtype=np.complex128)
    L = kron(I_ado, I_sys, format='csr')*0+0.j

    # System free liouviliian matrix
    H= sim.pot.H_mat
    I= np.eye(sim.pot.ns)
    Lsys0 = csr_matrix(-1.j*(np.kron(H,I) - np.kron(I,H.T)), dtype=np.complex128)
    L += kron(I_ado, Lsys0, format='csr')

    # Add on Xi
    Xi= sim.Xi
    if sim.params.LTCorr == 'same_for_each_ADO':
        L += kron(I_ado, Xi, format='csr')
    elif sim.params.LTCorr == 'different_for_each_ADO':
        print('not implemented yet')

    # Add on the off-diagonal coupling
    V= sim.pot.s_mat
    VL= kron(V,I, format='csr')
    VR= kron(I,V.T, format='csr')
    Vx= VL-VR
    for ki in range(K):
        Ck=sim.bath.C_ks[ki]
        Ak=A[ki]
        Adagk=Adag[ki]
        Lk_plus = -(1.j/np.sqrt(np.abs(Ck)))*(Ck*VL-Ck.conj()*VR)
        Lk_minus = -1.j*np.sqrt(np.abs(Ck))*Vx
        L += kron(Adagk, Lk_plus, format='csr')
        L += kron(Ak, Lk_minus, format='csr')

    # Add on the diagonal damping
    for ki in range(K):
        gamk=sim.bath.gam_ks[ki]
        AdagAk=AdagA[ki]
        L += -gamk*kron(AdagAk, I_sys, format='csr')


    # set the outputs
    sim.Liouvillian = L
    sim.params.Nados=Nados
    sim.params.Ntot=sim.params.Nados*sim.params.ns**2
    
