
import numpy as np
import sys


# PADE decomposition from the Hu et al. papaer

# f = e^[x/2] / (e^[x/2]-e^[-x/2])
#f-1/2 =  (e^[x/2]-  1/2(e^[x/2]-e^[-x/2])) / (e^[x/2]-e^[-x/2])
# = 1/2 (e^[x/2]+e^[-x/2]) / (e^[x/2]-e^[-x/2])
# = 1/2 coth(x/2)

def get_coeffs(mode,N,terminate=False):
    '''
    Get the coefficients for the Pade decomposition of the coth function
    mode = '[N/N]' or '[N-1/N]' for different types of Pade decomposition
    N = number of poles
    terminate = are we using a terminator?
    if so, we add lots of extra poles
    '''
    if terminate: N=20
    mode= str(mode.strip()) #clean up the mode string
    if mode == '[N/N]': #just do [N/N]
        M=2*N+1                     #(4c)
        def b(n): return 2*n+1      #(4b)
        def delta(i,j): return 1 if i == j else 0
        Lambda = np.zeros((2*N+1,2*N+1))
        Lambda_tilde = np.zeros((2*N,2*N))
        for n in range(2*N+1):
            for m in range(2*N+1):
                n_ = n+1
                m_= m+1
                if n<2*N and m <2*N:
                    Lambda_tilde[n,m] += delta(m_,n_+1)/np.sqrt(b(m_+1)*b(n_+1))
                    Lambda_tilde[n,m] += delta(m_,n_-1)/np.sqrt(b(m_+1)*b(n_+1))
                Lambda[n,m] += delta(m_,n_+1)/np.sqrt(b(m_)*b(n_))
                Lambda[n,m] += delta(m_,n_-1)/np.sqrt(b(m_)*b(n_))
        # Get the eigenvalues
        eigvals_Lambda,_=np.linalg.eigh(Lambda)
        eigvals_Lambda_tilde,_=np.linalg.eigh(Lambda_tilde)
        # Calculate the zetas and xis
        zeta = 2./np.flip(eigvals_Lambda_tilde)[0:N]
        xi =  2./np.flip(eigvals_Lambda)[0:N]
        # Calculate the residues
        eta=np.zeros(N)
        for j in range(N):
            num_prod = 1; denom_prod = 1
            for k in range(N):
                num_prod *= (zeta[k]**2 - xi[j]**2)
                if k!=j: denom_prod *= (xi[k]**2 - xi[j]**2)
            eta[j] = num_prod/denom_prod
        R_N = 1/(4*(N+1)*b(N+1))
        eta*= R_N/2
        return eta, xi, R_N, N #residues, poles, constant factor, number of poles

    elif mode=='[N-1/N]':
        M=2*N                         #(4c)
        def b(n): return 2*n+1      #(4b)
        def delta(i,j): return 1 if i == j else 0
        Lambda = np.zeros((2*N,2*N))
        Lambda_tilde = np.zeros((2*N-1,2*N-1))
        for n in range(2*N):
            for m in range(2*N):
                n_ = n+1
                m_= m+1
                if n<2*N-1 and m <2*N-1:
                    Lambda_tilde[n,m] += delta(m_,n_+1)/np.sqrt(b(m_+1)*b(n_+1))
                    Lambda_tilde[n,m] += delta(m_,n_-1)/np.sqrt(b(m_+1)*b(n_+1))
                Lambda[n,m] += delta(m_,n_+1)/np.sqrt(b(m_)*b(n_))
                Lambda[n,m] += delta(m_,n_-1)/np.sqrt(b(m_)*b(n_))
        # Get the eigenvalues
        eigvals_Lambda,_=np.linalg.eigh(Lambda)
        eigvals_Lambda_tilde,_=np.linalg.eigh(Lambda_tilde)
        # Calculate the zetas and xis
        zeta = 2./np.flip(eigvals_Lambda_tilde)[0:N-1]
        xi =  2./np.flip(eigvals_Lambda)[0:N]
        # Calculate the residues
        eta=np.zeros(N)
        for j in range(N):
            num_prod = 1; denom_prod = 1
            for k in range(N):
                if k<N-1: num_prod *= (zeta[k]**2 - xi[j]**2)
                if k!=j: denom_prod *= (xi[k]**2 - xi[j]**2)
            eta[j] = num_prod/denom_prod
        eta*=(N/2)*b(N+1)
        return eta, xi, 0, N #residues, poles, constant factor, number of poles
    elif mode=='[N+1/N]':
        sys.exit('not implemented yet')
    else: sys.exit('invalid mode')

    return


if __name__ == '__main__':
    N = 10
    import matplotlib.pyplot as plt
    import sys
    eta,xi,R_N = get_coeffs('[N-1/N]',N)
    # Plot the exact function vs the Pade approximation
    def xPhi_2Np1_x2(x):
        result= R_N*x
        for j in range(0,len(eta)):
            result+= (2*eta[j]*x/(x**2+xi[j]**2))
        return result

    def Phi_2Np1_x2(x):
        result= R_N
        for j in range(0,len(eta)):
            result+= (2*eta[j]/(x**2+xi[j]**2))
        return result


    def pade_NoN(x):
        return 1/x + 1/2 + xPhi_2Np1_x2(x)


    x = np.linspace(-300,300,10000)
    f_bose_exact= 1/(1-np.exp(-x))


    ax=plt
    # fig, (ax) = plt.figure()
    # ax.plot(x, f_bose_exact, label='exact results')
    # ax.plot(x, pade_NoN(x), label='approximate results')
    ax.plot(x, Phi_2Np1_x2(x), label='approximate results')
    ax.plot(x, (1/(2*np.tanh(x/2))-1/x)/x, label='exact results (with a coth)')
    # ax.ylim(-10,10)
    plt.legend()
    plt.show()

