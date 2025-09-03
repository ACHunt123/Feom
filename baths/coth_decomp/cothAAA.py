import os
import sys
import numpy as np
# Get the directory of the current script file
script_dir = os.path.dirname(os.path.abspath(__file__))
# go back one to get the bath location
script_dir = os.path.dirname(script_dir)


def get_coeffs(params, support=None, values=None,terminate=False):
    terminate=False # get the correct number of poles
    # Bathmode and settings
    minres_tol = 1e-6      # tolerance for the minimum abs value of a residue in the AAA decomposition

    ## Calculate the proposed extent of the support such that J(w) has decayed to 0
    if support is None or values is None: 
        Jw_min_tol = 1e-3      # tolerance for the maximum frequency of the grid for the AAA decomposition
        w_max=params.gam # start with the cuttoff frequency
        while params.J(w_max) > Jw_min_tol: w_max += 10 # find the maximum frequency where J(w) is still non-zero
        print(f'Calculating the support and values for the AAA decomposition with tolerance {Jw_min_tol}...')
        nw=int(w_max*10)
        if nw < 1000: nw = 1000  # ensure that the support is not too small
        support = np.linspace(-w_max,w_max,nw,dtype=np.complex128) # support for the AAA decomposition
        values = params.P(support)                     # values of the pole function at the support points
        print(f'Maximum frequency for the AAA decomposition: {w_max} (tolerance {Jw_min_tol})')
        ext = f'_nw{nw}_wmax{int(w_max)}.txt'                        # extension for the aaa files
    else:
        nw= len(support)  # number of support points
        print(f'Using provided support and values for the AAA decomposition.')
        ext = f'_quadrature_nw{nw}.txt'

    ### Use the AAA decomposition to get the coefficients
    mu_eff = params.mu if not terminate else 0  # number of poles for AAA decomposition (0 means as many as needed)
    folder = f'aaa_K{mu_eff}'                                   # folder to save the aaa files
    aaa_filename = f'aaa_data{ext}'                             # filename to save the aaa support and values data
    aaa_data_path = f'{folder}/{aaa_filename}'
    command= f"run_aaa_fromfile({mu_eff},'{os.getcwd()}/{folder}','{os.getcwd()}/{aaa_data_path}','{ext}',{str(terminate).lower()})" # the command to run the AAA decomposition in MATLAB

    if not os.path.exists(folder): os.makedirs(folder)
    if not os.path.exists(aaa_data_path):  # save the support and values to a file if it does not exist
        data = np.column_stack((support.real, values.real, values.imag))  
        print(f'Saving support and values to {aaa_data_path} ...')
        np.savetxt(f'{aaa_data_path}', data, header='Support Re[Values] Im[Values]', comments='')         # save the support and values to be read by matlab
    ### Run the AAA decomposition in MATLAB
    if(0):                                      # run the MATLAB script using the system command
        print('Running AAA decomposition in MATLAB...')
        os.system(f"matlab -batch 'cd {script_dir}/aaa;{command}' > /dev/null 2>&1")    # run the matlab script to get the AAA coefficients
        print('AAA decomposition complete, loading results...')
    else:                                       # print out the command to run the MATLAB script if it has not already been run
        if not os.path.exists(f'{folder}/pol_real{ext}'):
            print('run the following command in MATLAB to get the AAA coefficients:')
            print(f"\n{command}")
            ### append the command to a file for later use
            with open(f'{script_dir}/aaa/commands_to_run.m', 'a') as f:
                if params.mu>0: f.write(f"{command}\n")
            print(f"\n")
            sys.exit()
        else:
            print('AAA decomposition already done, loading results from files...')

    ### Load the aaa results
    repoles = np.loadtxt(f'{folder}/pol_real{ext}')
    impoles = np.loadtxt(f'{folder}/pol_imag{ext}')
    params.poles = repoles + 1.j * impoles
    reres = np.loadtxt(f'{folder}/res_real{ext}')
    imres = np.loadtxt(f'{folder}/res_imag{ext}')
    konstant = np.loadtxt(f'{folder}/k{ext}')  # load the constant shift k
    params.res = reres + 1.j * imres
    params.res_original = params.res.copy()         # save the original residues for later use
    params.poles_original = params.poles.copy()     # save the original poles for later use
    ### Clean up the poles and residues
    mask= np.abs(params.res)> minres_tol
    params.res = params.res[mask]      # remove any tiny residues
    params.poles = params.poles[mask]  # remove the corresponding poles
    params.res = np.imag(params.res)*1.j            # remove the real parts, as by symmetry they should be zero
    params.poles = np.imag(params.poles)*1.j
    ### Calulate the real coefficients w_i and gamma_i from conjugate pairs of poles and residues
    upper_poles= []; upper_res = []
    for k in range(len(params.poles)):
        if np.imag(params.poles[k]) > 0:
            upper_poles.append(params.poles[k])
            upper_res.append(params.res[k])
    upper_poles = np.array(upper_poles,dtype=np.complex128) ; upper_res = np.array(upper_res,dtype=np.complex128)
    w_i = np.imag(upper_poles)                             # these are the new prequencies
    gam_i = -2*np.imag(upper_poles)*np.imag(upper_res)     # these are the new gammas
    return gam_i, w_i, konstant, len(gam_i)  # return the gammas, frequencies, constant shift k, and number of exponentials