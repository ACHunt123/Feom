#!/usr/bin/env python3
import numpy as np
import sys
import matplotlib.pyplot as plt
# This file contains the coefficients for the poles of the bath correlation function

data=[[1,-0.007050248,0.002482056,188.7752364,-164.9169644],
[2,0.102762949,0.116652539,114.5230822,-130.046728],
[3,0.581303135,-0.644477046,69.37948966,-100.2912245],
[4,-1.455575575,-1.650266045,40.55914002,-74.97480723],
[5,-2.988444895,1.187910091,22.52508418,-54.17566836],
[6,-0.622081237,3.115931332,11.81278425,-37.82551602],
[7,1.30388447,2.013121889,5.852688004,-25.57526889],
[8,1.382279145,0.625686611,2.752595646,-16.83290689],
[9,0.8611309,0.016022,1.241895418,-10.8518732],
[10,0.443366368,-0.118339778,0.548086947,-6.885813753],
[11,0.212002544,-0.101557778,0.243348806,-4.31568345],
[12,0.09893136,-0.063028266,0.113960488,-2.680510104],
[13,0.046158834,-0.034415289,0.059450913,-1.654822486],
[14,0.021773639,-0.017637193,0.033910866,-1.017022414],
[15,0.010348793,-0.008760985,0.019524006,-0.62200475],
[16,0.004904562,-0.004276533,0.010782038,-0.378232958],
[17,0.002297321,-0.002058266,0.005782001,-0.228876357],
[18,0.001065577,-0.00097944,0.003293429,-0.138297459],
[19,0.00049681,-0.000464088,0.002202629,-0.083659132],
[20,0.000235641,-0.000219394,0.001675626,-0.050533949],
[21,0.000112977,-0.000102543,0.001268055,-0.03029634],
[22,5.37E-05,-4.70E-05,0.000861682,-0.017943686],
[23,2.49E-05,-2.12E-05,0.000496313,-0.010482398],
[24,1.12E-05,-9.49E-06,0.000227853,-0.006041449],
[25,4.85E-06,-4.23E-06,6.86E-05,-0.003437852],
[26,2.05E-06,-1.89E-06,-9.00E-06,-0.001928643],
[27,8.36E-07,-8.31E-07,-3.82E-05,-0.001060056],
[28,3.24E-07,-3.56E-07,-4.17E-05,-0.000568788],
[29,1.21E-07,-1.49E-07,-3.41E-05,-0.000297603],
[30,4.48E-08,-6.28E-08,-2.43E-05,-0.000147216],
[31,1.62E-08,-2.71E-08,-1.54E-05,-5.76E-05]]
# Convert data to numpy array
data = np.array(data)
# Extract columns
K = data.shape[0]  # Number of poles
K = 30  # Number of poles
n = data[:K, 0].astype(int)  # Pole number
eta = data[:K, 1] + 1.j* data[:K,2]               # Coefficient residue
xi = data[:K, 3] + 1.j* data[:K,4]             # Pole frequency

# define S(omega)
def S(omega):
    one = np.ones_like(eta)
    return np.sum(eta / (omega*one - xi)) + np.sum(eta.conj() / (omega*one - xi.conj()))

if(0): #plot the poles and the function S(ω)
    # Frequency range
    w = np.linspace(-300, 300, 10000)  # Frequency range
    S = np.array([S(omega) for omega in w])  # Calculate S(ω) for each frequency
    # plot the poles
    # make two axes, one for the poles and one for the function S(ω)
    fig, (pole_ax, func_ax) = plt.subplots(1, 2, figsize=(6.5, 3.2))  # width, height in inches

    # Set font size to match LaTeX \normalsize (about 10pt)
    plt.rcParams.update({'font.size': 10})
    pole_ax.set_title(r'Poles of $S(\omega)$')
    pole_ax.set_xlabel('Real Part')
    pole_ax.set_ylabel('Imaginary Part')
    # pole_ax.grid(True)
    pole_ax.scatter(xi.real, xi.imag, color='blue', label='Poles',s=2)
    pole_ax.scatter(xi.real, -xi.imag, color='blue',s=2)  # Mirror image for symmetry

    #add an inset fig with a zoomed-in view of the poles
    inset_ax = pole_ax.inset_axes([0.525, 0.3, 0.4, 0.4])  # [x, y, width, height]
    inset_ax.scatter(xi.real, xi.imag, color='blue',s=2)  # Poles
    inset_ax.scatter(xi.real, -xi.imag, color='blue',s=2)  # Mirror image for symmetry
    inset_ax.set_ylim(-0.0005, 0.0005)  # Adjust limits for zoom
    inset_ax.set_xlim(-0.0005, 0.0005)
    inset_ax.axhline(0, color='black', lw=0.5, ls='--')
    inset_ax.set_title('Zoomed-in',fontsize=8)
    # make the tick labels smaller
    for label in inset_ax.get_xticklabels() + inset_ax.get_yticklabels():
        label.set_fontsize(7)

    # plot with dashed lines
    func_ax.plot(w, S.real, label=r'$S(\omega)$', color='red', linestyle='-')
    func_ax.plot(w, S+np.flip(S), label=r'$S^{+}(\omega)$', color='blue', linestyle='-.')
    func_ax.plot(w, S-np.flip(S), label=r'$S^{-}(\omega)$', color='orange', linestyle=':')
    func_ax.plot(w, (S+np.flip(S))/(S-np.flip(S)), label=r'$\frac{S^{+}(\omega)}{S^{-}(\omega)}\equiv A(\omega)$', color='purple', linestyle='-.')
    # func_ax.plot(w, S.imag, label='Im[S(ω)]', color='green')
    # plt.axhline(0, color='black', lw=0.5, ls='--')
    # plt.axvline(0, color='black', lw=0.5, ls='--')
    func_ax.set_xlabel('Frequency (ω)')
    func_ax.set_ylabel('S(ω)')
    #make legend fontsize smaller

    func_ax.legend(fontsize=8, loc='lower right')
    func_ax.set_title(r'$S(\omega)$ and derived functions')
    # plt.ylim(-5, 2)
    fig.tight_layout()
    # plt.title(r'$S', fontsize=14)
    plt.savefig("FreePoles_Sw_fromFig1.pdf", format="pdf")  # saved in the current working directory
    plt.show()

# Calculate fac*J(ω) from S(ω) 
def J(omega):
    return S(omega) - S(-omega)
# note that there are extra factors here = 4\pi\hbar,  but who cares as we just want to find the zeros of J(ω)

# Output J(ω) for a range of frequencies and save for a MATLAB reading
if(0):
    w = np.linspace(-300, 300, 1000000)  # Frequency range
    J_values = np.array([J(omega) for omega in w])  # Calculate J(ω) for each frequency
    # Plot J(ω)
    plt.figure(figsize=(6.5, 3.2))  # width, height
    plt.plot(w, J_values.real, label=r'Re[J(ω)]', color='red', linestyle='-')
    plt.plot(w, J_values.imag, label=r'Im[J(ω)]', color='green', linestyle='--')
    plt.axhline(0, color='black', lw=0.5, ls='--')
    plt.axvline(0, color='black', lw=0.5, ls='--')
    plt.xlabel('Frequency (ω)')
    plt.ylabel(r'J(ω)')
    plt.title(r'J(ω) from Free Poles')
    plt.legend(fontsize=8, loc='upper right')
    plt.tight_layout()
    plt.show()
    # Save to file that can be read by matlab, and shift them so they are 0 at ends
    J_end = J_values[-1]  # Get the last value of J(ω)
    J_values -= J_end  # Shift to make the first value zero
    data= np.column_stack((w, J_values.real, J_values.imag))
    np.savetxt("FreePoles_Jw_fromFig1.txt", data, header="Frequency (ω) Re[J(ω)] Im[J(ω)]", delimiter='\t')
    print("Saved J(ω) values to FreePoles_Jw_fromFig1.txt")

# Read the poles and residues from the file from the AAA output of the previous step
if(0):
    ### Load the aaa files
    folder='.files'
    repoles = np.loadtxt(f'{folder}/pol_real.txt')
    impoles = np.loadtxt(f'{folder}/pol_imag.txt')
    poles = repoles + 1.j * impoles
    reres = np.loadtxt(f'{folder}/res_real.txt')
    imres = np.loadtxt(f'{folder}/res_imag.txt')
    res = reres + 1.j * imres


    # Setup the pole plot
    fig, (pole_ax, func_ax) = plt.subplots(1, 2, figsize=(12, 6))
    pole_ax.set_title('Poles')
    pole_ax.set_xlabel('Real Part')
    pole_ax.set_ylabel('Imaginary Part')
    pole_ax.grid(True)
    # plot the poles of the aaa output
    pole_ax.plot(repoles, impoles, 'o', label=f'AAA, N={len(poles)}', color='blue') 
    #plot the original poles
    pole_ax.plot(xi.real, xi.imag, 'x', label='Free Poles', color='orange')
    pole_ax.plot(-xi.real, -xi.imag, 'x', label='Free Poles', color='orange')
    pole_ax.plot(xi.real, -xi.imag, 'x', label='Free Poles', color='orange')
    pole_ax.plot(-xi.real, xi.imag, 'x', color='orange')  # complex conjugates

    # Frequency range
    plt.show()

# Calculate the coeficients for the Polynomials that correspond to the  simple pole functions
# Calculate the polynomial coefficients
pols=np.concatenate([xi.copy(), xi.conj().copy()])       # Include complex conjugates
reds=np.concatenate([eta.copy(), eta.conj().copy()])    # Include complex conjugates
pols = np.array(pols, dtype=complex)  # Ensure poles are complex
reds = np.array(reds, dtype=complex)  # Ensure residues are complex

def fac_prods_old(z):
    ''' Calculates the coefficients of w for the polynomial
    P(x) = prod_{i=1}^{n} (x + z_i)
    takes in the list of z
    and returns coefficients of the polynomial in the form:
    P(x) = c_0 + c_1*x + c_2*x^2 + ... + c_n*x^n
    where n is the number of elements in z
    '''
    v = np.zeros_like(z, dtype=int) # binary vector, which gives whether x or w is being multiplied in each bracket
    coeffs = np.zeros(len(z)+1, dtype=np.dtype(z[0]))  # Coefficients for the polynomial
    def binary_increment(v):
        ''' Inrements the binary vector by 1'''
        n_10 = np.sum(v*2**np.arange(len(v)-1, -1, -1))  # Convert binary vector to decimal
        n_10 += 1
        if n_10%100000==0:  # Print progress every 100 increments
            print(f"Progress: {n_10}/{2**len(z)}")
        v = np.array([int(i) for i in np.binary_repr(n_10, width=len(z))])
        return v
    go= True                # Continue until all combinations are processed
    while go==True:         #   continue untill all elements are 1 (ie we have gone through all combinations)
        order = np.sum(v)   # exponent of the polynomial (as there are)
        prod = 1
        for i in range(len(z)):
            if not v[i]: # If this z_i is selected for this bracket (not x)
                prod *= z[i]
        coeffs[order] += prod 
        # print(f"Current coeffs: {coeffs}, v: {v} prod {prod}")
        if np.all(v==1):
            go = False
        v = binary_increment(v)  # Increment the binary vector
    return coeffs

def product_of(coeffs_1, coeffs_2, len_out=None):
    ''' Multiplies two polynomials given their coefficients
    coeffs_1: Coefficients of the first polynomial
    coeffs_2: Coefficients of the second polynomial
    Returns the coefficients of the product polynomial
    All in form P(x) = c_0 + c_1*x + c_2*x^2 + ... + c_n*x^n
    '''
    n = len(coeffs_1) + len(coeffs_2) - 1  if len_out == None else len_out# Degree of the product polynomial
    product_coeffs = np.zeros(n, dtype=np.dtype(coeffs_1[0]))  # Initialize product coefficients
    for i in range(len(coeffs_1)):
        for j in range(len(coeffs_2)):
            if i+j==n: break
            product_coeffs[i + j] += coeffs_1[i] * coeffs_2[j]
    return product_coeffs

def fac_prods(z):
    ''' Calculates the coefficients of w for the polynomial
    P(x) = prod_{i=1}^{n} (x + z_i)
    takes in the list of z
    and returns coefficients of the polynomial in the form:
    P(x) = c_0 + c_1*x + c_2*x^2 + ... + c_n*x^n
    where n is the number of elements in z
    '''
    coeffs = np.zeros(len(z)+1, dtype=np.dtype(z[0]))  # Coefficients for the polynomial
    print(f"Calculating coefficients for {len(z)} poles")
    coeffs[0] = 1  # start with 1 + 0x + 0x^2 + ... + 0x^n    
    for z_i in z: 
        term = np.zeros_like(coeffs)
        term[0] = z_i
        term[1] = 1  # term = z_i + x
        coeffs = product_of(coeffs, term, len_out=len(coeffs))  # Multiply the current coefficients with the new term
    return coeffs


# if(0):  # Test the fac_prods function
#     z= np.array([1,266+2.j,1])
#     def fx(x): return np.prod(x + z)
#     cs = fac_prods(z)
#     Fx_expanded = np.poly1d(cs[::-1])  # Reverse the coefficients for np.poly1d
#     # plot the polynomial
#     x = np.linspace(-3, 3, 1000)
#     y = Fx_expanded(x)
#     plt.figure(figsize=(6.5, 3.2))  # width, height
#     plt.plot(x,  Fx_expanded(x), label='Polynomial from expansion', color='blue')
#     plt.plot(x,  [fx(xi) for xi in x], label='Polynomial ', color='red', linestyle='--')
#     plt.xlabel('x')
#     plt.ylabel('P(x)')
#     plt.title('Polynomial from Free Poles')
#     plt.show()
#     sys.exit()

def S_coefs(poles, residues):
    ''' Calculate the coefficients of the numerator and denominator of S(w) from the poles and residues
    F(x) =  sum_{all i=1}^n eta_i prod_{all i'!=i}^{n} (x - /x_i')
    where eta_i are the poles and xi_i are the residues, and the sum over `all i' means sum over conjugates also.

    Input only one of the conjugation of the poles and residues, we will calculate in steps:

    1. Calculate the coefficients of the numerator and denominators for
    sum_{i=1}^n eta_i/(x - xi_i) = P(x)/Q(x)
    where P(x) = sum_{i=1}^n eta_i prod_{i'!=i}^{n} (x - /xi_i')
    and Q(x) = prod_{i=1}^n (x - xi_i)
    and the sums are only over a single set of poles and residues.
    F(x) = P(x)/Q(x) + P*(x)/Q*(x)

    2. Calculate the numerator coefficients for S(w)
    = P(x)Q*(x) + P*(x)Q(x)

    3. Calculate the denominator coefficients for S(w) 
    = Q(x)Q*(x) 
    '''
    # 1.
    assert poles.size == residues.size, "\nPoles and residues must have the same size"
    assert len(poles) > 3, "\nPoles array must be longer than 3 for the slicing processes to work"
    P_coeffs = np.zeros(len(poles) + 1, dtype=np.dtype(poles[0]))  # Coefficients for P(x)
    Q_coeffs = np.zeros(len(poles) + 1, dtype=np.dtype(poles[0]))  # Coefficients for Q(x)
    for i,eta_i in enumerate(residues):
        print(i/len(residues), end='\r')  # Print progress
        xi_i_prime = np.delete(poles.copy(), i)  # Remove the i-th pole
        print(eta_i)
        print(xi_i_prime)
        print(fac_prods(-xi_i_prime))
        P_coeffs[:-1] += eta_i * fac_prods(-xi_i_prime)
        if i==0:    #save the expansion of the first bracket for the next step
            Q_coeffs_no_xi0 = P_coeffs/eta_i # Coefficients of the denominator without the first pole
            xi_0 = poles[0]  # The first pole
    # Add the first pole back to the denominator coefficients
    Q_coeffs[1:] = Q_coeffs_no_xi0[:-1]
    Q_coeffs = Q_coeffs - xi_0*Q_coeffs_no_xi0
    # 2.
    Sw_numerator_coeffs = product_of(P_coeffs, Q_coeffs.conj()) + product_of(P_coeffs.conj(), Q_coeffs)  # Coefficients for the numerator of S(w)
    # 3.
    Sw_denominator_coeffs = product_of(Q_coeffs, Q_coeffs.conj())  # Outer product for P(x)Q*(x) 

    # Return the coefficients of the numerator and denominator
    return Sw_numerator_coeffs, Sw_denominator_coeffs

if(1): #test the S_coefs function
    # test poles
    pols = np.array([1+2.j, 266 + 2.j, 1+3j,2.j])  # Example poles
    reds = np.array([1, 1, 1,4.j])  # Example residues

    pols = xi.copy()  # Use the poles from the data
    reds = eta.copy()  # Use the residues from the data
    # Function to evaluate S(w)
    def Sw_exact(xi, eta, omega):
        return np.sum(eta / (omega - xi)) + np.sum(eta.conj() / (omega - xi.conj()))
    # Calculate the coefficients of the numerator and denominator
    Numerator_coeffs_list, Denominator_coeffs_list = S_coefs(pols, reds)
    
    x = np.linspace(-300, 300, 1000)
    # Calculate the polynomial from the coefficients
    y_exact = [Sw_exact(pols, reds, xi) for xi in x]  # Calculate S(ω) for each frequency
    y_poly = np.polyval(Numerator_coeffs_list[::-1], x) / np.polyval(Denominator_coeffs_list[::-1], x)
    # Plot the results
    plt.figure(figsize=(6.5, 3.2))  # width, height
    plt.plot(x, np.real(y_exact), label='Exact S(ω)', color='red')
    plt.plot(x, np.imag(y_exact), label='Exact S(ω)', color='green')
    plt.plot(x, np.real(y_poly), label='Polynomial from Free Poles', color='blue', linestyle='--')
    plt.plot(x, np.imag(y_poly), label='Polynomial from Free Poles', color='orange', linestyle='--')
    plt.show()
    sys.exit()
# plot the polynomial
x = np.linspace(-300, 300, 1000)
y = np.polyval(Numerator_coeffs_list, x)/np.polyval(Denominator_coeffs_list, x)
plt.figure(figsize=(6.5, 3.2))  # width, height
plt.plot(x, np.real(y), label='Polynomial from Free Poles', color='blue')
plt.plot(x, np.imag(y), label='Polynomial from Free Poles', color='green')
plt.plot(x,np.real([(Swtemp(pols,reds,x_i))for x_i in x]), label='S(ω) from Free Poles', color='red', linestyle='--')
plt.ylim(-10, 10)  # Set y-limits for better visibility
plt.figure(figsize=(6.5, 3.2))  # width, height
plt.plot(np.real(pols), np.imag(pols), 'o', label='Poles', color='blue')
plt.show()


