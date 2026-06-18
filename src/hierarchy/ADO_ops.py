import numpy as np
from scipy.sparse import eye, kron
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
    rho0[ns**2:2*ns**2]=(np.identity(ns)).flatten()

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
    print(rho0[:ns**2].reshape((ns,ns)))
    for k in range(K+3):
        # rho0 =Adagks.dot(rho0)
        # rho0 =AdagAks.dot(rho0)
        rho0 =Aks.dot(rho0)
        print(f'Ak ^{k+1}rho000...')
        print(rho0[:ns**2].reshape((ns,ns)))
        # print(rho0[ns**2:2*ns**2].reshape((ns,ns)))
        
        print('the whole matrix is zero?',np.allclose(rho0 , 0))

    # this shows that the Ak will bring down the rho_n ->rho_n-1
    # meaning that rho_n will interact with rho_n+1