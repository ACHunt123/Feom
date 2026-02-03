import os
import sys
import numpy as np
import Feom.baths.aaa.pyA4 as pyA4
# Get the directory of the current script file
script_dir = os.path.dirname(os.path.abspath(__file__))
# go back one to get the bath location
script_dir = os.path.dirname(script_dir)


def get_coeffs(params, support, values, max_accuracy=False,ext_fname=''):
    max_accuracy=False # get the correct number of poles
    # Bathmode and settings
    minres_tol = 1e-6               # tolerance for the minimum abs value of a residue in the AAA decomposition
    ext = getattr(params, 'support_param_str', '')  + ext_fname + '.txt' # extension for the aaa files

    ### Use the AAA decomposition to get the coefficients
    mu_eff = params.mu if not max_accuracy else 0  # number of poles for AAA decomposition (0 means as many as needed)
    folder = f'aaa_K{mu_eff}{ext_fname}'                                   # folder to save the aaa files
    aaa_filename = f'aaa_data_{ext}'                             # filename to save the aaa support and values data
    aaa_data_path = f'{folder}/{aaa_filename}'
    command= f"run_aaa_fromfile({mu_eff},'{os.getcwd()}/{folder}','{os.getcwd()}/{aaa_data_path}','_{ext}',{str(max_accuracy).lower()})" # the command to run the AAA decomposition in MATLAB

    if not os.path.exists(folder): os.makedirs(folder)
    if not os.path.exists(aaa_data_path):  # save the support and values to a file if it does not exist
        data = np.column_stack((support.real, values.real, values.imag))  
        print(f'Saving support and values to {aaa_data_path} ...')
        np.savetxt(f'{aaa_data_path}', data, header='Support Re[Values] Im[Values]', comments='')         # save the support and values to be read by matlab
        
    ### Run the AAA decomposition in python (if it hasnt already been done)
    if not os.path.exists(f'{folder}/pol_real_{ext}'): # check if the outfiles exist already
        pyA4.run_aaa_fromfile(mu_eff,f'{os.getcwd()}/{folder}',f'{os.getcwd()}/{aaa_data_path}',f'_{ext}',max_accuracy)
        print('AAA decomposition complete, loading results...')
    else:
        print('AAA decomposition already done, loading results from files...')
    ### Load the aaa results
    repoles = np.loadtxt(f'{folder}/pol_real_{ext}')
    impoles = np.loadtxt(f'{folder}/pol_imag_{ext}')
    params.poles = repoles + 1.j * impoles
    reres = np.loadtxt(f'{folder}/res_real_{ext}')
    imres = np.loadtxt(f'{folder}/res_imag_{ext}')
    konstant = np.loadtxt(f'{folder}/k_{ext}')  # load the constant shift k
    params.res = reres + 1.j * imres
    params.res_original = params.res.copy()         # save the original residues for later use
    params.poles_original = params.poles.copy()     # save the original poles for later use
    
    gam_i = np.atleast_1d(np.loadtxt(f'{folder}/gam_i_{ext}'))
    w_i = np.atleast_1d(np.loadtxt(f'{folder}/w_i_{ext}'))
    konstant = np.loadtxt(f'{folder}/k_{ext}')  # load the constant
    return gam_i, w_i, konstant, len(gam_i)   # return the gammas, frequencies, constant shift k, and number of exponentials


def markovian_pole_trunc(params, max_freq=2.5e2):
    ''' Remove any poles with frequencies above max_freq  and add their contribution to the constant shift k'''
    propagate = np.abs(params.gam_ks) <= max_freq
    markovian = np.abs(params.gam_ks) > max_freq

    mark_gam_ks = params.gam_ks[markovian]
    mark_C_ks = params.C_ks[markovian]

    # assert that the coefficients are real
    assert np.all(np.isreal(mark_C_ks)), 'Markovian coefficients are not real'
    assert np.all(np.isreal(mark_gam_ks)), 'Markovian gammas are not real'

    # calculate the extra constant shift k from the markovian poles
    extra_k = (params.beta/(2*params.eta*params.gam)) * np.real(np.sum(mark_C_ks/mark_gam_ks))
    
    # update the constant shift k and the frequencies and coefficients
    params.k += extra_k
    params.gam_ks = params.gam_ks[propagate]
    params.C_ks = params.C_ks[propagate]
    # update the number of exponentials
    params.N_exp = len(params.gam_ks)
    params.K = params.N_exp - params.N_nonmats
    
    print(f'Removed {np.sum(markovian)} high frequency poles above from the bath, with original number of poles {len(params.gam_ks)+len(mark_gam_ks)}')
    return