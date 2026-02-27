import numpy as np
from scipy.sparse import csr_matrix, eye, kron
from scipy.sparse import coo_matrix
from types import SimpleNamespace


def generate_ado_raising_lowering_ops(K,L):
    # Initialize lists of rows and columns for the lowering matrix
    A_rows = [[] for _ in range(K)]; A_cols = [[] for _ in range(K)]; A_data = [[] for _ in range(K)]
    indices_dict = {tuple('0' * K): 0}
    old_indices_dict = {tuple('0' * K): 0}

    for tier in range(L):
        new_indices_dict = {}
        for old_tuple, old_id in old_indices_dict.items():
            for k in range(K):
                # Extract the CURRENT excitation number for this specific mode
                current_n_k = int(old_tuple[k])
                new_tuple = list(old_tuple)
                new_tuple[k] = str(current_n_k + 1)
                new_tuple = tuple(new_tuple)
                
                # Register if it's new
                if new_tuple not in indices_dict:
                    new_id = len(indices_dict)
                    indices_dict[new_tuple] = new_id
                    new_indices_dict[new_tuple] = new_id
                else:
                    new_id = indices_dict[new_tuple] # Fetch the existing ID!

                # The transition value depends strictly on the mode's excitation number
                val = np.sqrt(current_n_k + 1)
                # --- A_k (Lowering operator): Takes new_id (col) to old_id (row) ---
                A_rows[k].append(old_id)
                A_cols[k].append(new_id)
                A_data[k].append(val)
                # This way, for the example of one dim
                # <n| A \sum_n' \rho_n' |n'> 
                #  = \sum_n' \rho_n' <n| A |n'> 
                #  = \sum_n' \rho_n' delta_{n n'-1} 
                #  = rho_n+1 
        old_indices_dict = new_indices_dict.copy()

    # Compile into perfectly square K distinct sparse matrices
    Nados = len(indices_dict)
    A = [coo_matrix((A_data[k], (A_rows[k], A_cols[k])), shape=(Nados, Nados)).tocsr() for k in range(K)]
    Adag = [np.transpose(np.conj(A[k]))  for k in range(K)]
    AdagA = [Adag[k] @ A[k] for k in range(K)]


    # Package and send out the operators
    ADO_ops = SimpleNamespace(A=A,Adag=Adag,AdagA=AdagA)
    return ADO_ops,indices_dict


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
    

if __name__=='__main__':
    # test for the A and the A daggers
    ns = 2
    K=4
    L=3

    # get the superoperators
    I_sys = eye(ns**2, format='csr', dtype=np.complex128)
    ADO_ops,indices_dict = generate_ado_raising_lowering_ops(K,L)
    Nados = len(indices_dict)

    
    I_ado = eye(Nados, format='csr', dtype=np.complex128)

    A = ADO_ops.A
    Adag = ADO_ops.Adag
    AdagA = ADO_ops.AdagA
    Nados = len(indices_dict)

    rho0 = np.zeros(Nados*ns**2,dtype=complex)

    #set the 0th ado to nonzeros
    print('setting one of the rho100... to the identity and the others to 0')
    rho0[ns**2:2*ns**2]=(np.identity(ns)).flatten(order='F')

    # calculate the lowering and raising superoperator (summed over k)
    Aks = kron(I_ado, I_sys, format='csr')*0+0.j
    Adagks = kron(I_ado, I_sys, format='csr')*0+0.j
    AdagAks = kron(I_ado, I_sys, format='csr')*0+0.j
    for ki in range(K):
        Adagk=Adag[ki]
        AdagAk=AdagA[ki]
        Ak=A[ki]
        Aks += kron(Ak, I_sys, format='csr')
        Adagks += kron(Adagk, I_sys, format='csr')
        AdagAks += kron(AdagAk, I_sys, format='csr')

    print('original rho000...')
    print(rho0[:ns**2].reshape((ns,ns),order='F'))
    for k in range(K+3):
        # rho0 =Adagks.dot(rho0)
        # rho0 =AdagAks.dot(rho0)
        rho0 =Aks.dot(rho0)
        print(f'Ak ^{k+1}rho000...')
        print(rho0[:ns**2].reshape((ns,ns),order='F'))
        # print(rho0[ns**2:2*ns**2].reshape((ns,ns),order='F'))
        
        print('the whole matrix is zero?',np.allclose(rho0 , 0))

    # this shows that the Ak will bring down the rho_n ->rho_n-1
    # meaning that rho_n will interact with rho_n+1