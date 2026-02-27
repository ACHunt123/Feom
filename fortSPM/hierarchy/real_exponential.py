import numpy as np
from scipy.sparse import csr_matrix, eye, kron
from scipy.sparse import coo_matrix


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
                    
        old_indices_dict = new_indices_dict.copy()

    # Compile into perfectly square K distinct sparse matrices
    Nados = len(indices_dict)
    A = [coo_matrix((A_data[k], (A_rows[k], A_cols[k])), shape=(Nados, Nados)).tocsr() for k in range(K)]
    Adag = [np.transpose(np.conj(A[k]))  for k in range(K)]
    AdagA = [Adag[k] @ A[k] for k in range(K)]


    if(0):
        print(f"Total ADOs generated: {Nados}")
        # print(f"indices dict: {indices_dict}\n")

        for k in range(K):
            print(f"--- Adag Matrix for Mode k={k} ---")
            print(np.round(Adag[k].toarray(), 3))
            print(np.round(AdagA[k].toarray(), 3))
            print()
        exit()


    return A,Adag,AdagA,indices_dict


def generate_liouvillian(sim):

    K = len(sim.bath.C_ks)
    L = sim.params.L

    # Raising/Lowering and number operators in the ADO space
    A, Adag, AdagA, indices_dict = generate_ado_raising_lowering_ops(K,L)
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

        L += kron(Ak, Lk_plus, format='csr')
        L += kron(Adagk, Lk_minus, format='csr')

    # Add on the diagonal damping
    for ki in range(K):
        gamk=sim.bath.gam_ks[ki]
        AdagAk=AdagA[ki]
        L += -gamk*kron(AdagAk, I_sys, format='csr')


    # set the outputs
    sim.Liouvillian = L
    sim.params.Nados=Nados
    sim.params.Ntot=sim.params.Nados*sim.params.ns**2
