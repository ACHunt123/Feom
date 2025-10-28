% File: TFaySB.m
%%% Wrapper function for the spin boson model for Tom Fay's code
% to run in terminal, use:
% matlab -batch "addpath('~/software/phd/Feom/benchmark'); TFaySB(###)"
% where ### is the output folder path
%
function TFaySB(outfolder)

if nargin < 1
    outfolder = '';  % default
end
%%%

%%% path to the HEOMLAB library
addpath(genpath('/home/ach221/software/phd/HEOMLAB/heom-lab/functions'));
%%%

%%% Parameters for the problem
% system hamiltonian parameters (harmonic oscillator)
omega=1;
mass=1740;
% debye bath parameters
beta = 0.01 ;
lambda_D = 500 ; % should be half of mine
omega_D = 5.0 ; % should be the same as mine

% dynamics information - parameters for the Short-Iterative Arnoldi
% Integrator
dt = 1e-2 ;
n_steps = 1000 ;
krylov_dim = 8 ;
krylov_tol = 1e-8 ;
% parmeters for heirarchy truncation using L/M truncation
M_max = 2; 

% get the basis and operators for the harmonic oscillator
ns=10;
dx=0.01;
xmin=-5;
xmax=5;
hbar=1;
nx = floor((xmax - xmin) / dx);


[H0, pos_mat, x_arr, psi_ns, E_ns] = ho_matrices(mass, omega, hbar, xmin, dx, nx, ns);


% matrices of system observable operators to be returned, 
O_sys = {pos_mat} ;

% initial state of the system
delEs = E_ns - E_ns(1);         % energy differences from ground state
rho_s = diag(exp(-beta * delEs));   % Boltzmann weights
Zs = trace(rho_s);              % partition function
rho_0_sys = (rho_s * pos_mat) / Zs;

% renormalization potential (x2 as the thing is half mine, but the term added to the hamiltonian should be the same)
Hren=lambda_D.*pos_mat*pos_mat;

% two objects are supplied to the HEOM dynamics function:
% "full_system" specifies the full Hamiltonian (system + bath) and the
% temperature.
% "heom_dynamics" specifies the HEOM truncation, integrator for the
% dynamics and the total propagation time and observables to be calculated.

% the full_system object contains all information about the Hamiltonian of
% the full open quantum system
full_system = struct ;
% H_sys contains the system Hamiltonian
full_system.H_sys = H0+Hren;
full_system.beta = beta ;
% a struct that contains information about the HEOM dynamics
heom_dynamics = struct ;
% integrator information, currently only the short iterative arnoldi is
% implemented
heom_dynamics.integrator = struct ;
heom_dynamics.integrator.method = "SIA" ;
heom_dynamics.integrator.dt = dt ;
heom_dynamics.integrator.n_steps = n_steps ;
heom_dynamics.integrator.krylov_dim = krylov_dim ;
heom_dynamics.integrator.krylov_tol = krylov_tol ;

% hierarchy trunction information
heom_dynamics.heom_truncation = struct ;
heom_dynamics.heom_truncation.truncation_method = "depth cut-off" ;
heom_dynamics.heom_truncation.M_max = M_max ;
% heom_dynamics.heom_truncation.truncation_method = "frequency cut-off" ;
% heom_dynamics.heom_truncation.Gamma_cut = Gamma_cut ;


% what system observables should be returned
heom_dynamics.observables = struct ;
heom_dynamics.observables.system = O_sys ;

% set the initial condition
heom_dynamics.rho_0_sys = rho_0_sys ;

%% run the dynamics LOOPING OVER DIFFERENT PARAMETER SETS

L_max_list={2};
full_system_bathslist={struct("V",pos_mat,"spectral_density","debye","omega_D",omega_D,"lambda_D",lambda_D)};
termination_list={"low temp correction"};
namelist={'ITlowtemp'};




% extra terms for the NZ2 terminator 
heom_dynamics.heom_truncation.termination_k_max = 500 ; % max number of mats terms for NZ2
% heom_dynamics.heom_truncation.diagonal_only_term = false ;



% loop over the different L_max values
for i = 1:1
    L_max = L_max_list{i};
    heom_dynamics.heom_truncation.L_max = L_max;
    full_system.baths = {full_system_bathslist{i}} ;
    heom_dynamics.heom_truncation.heom_termination = termination_list{i};
    name=namelist{i};
  
    [O_t,t] = runHEOMDynamics(full_system,heom_dynamics) ;


    %%% save results to a text file
    %make the filename 
    % Create a descriptive filename
    filename = sprintf('TFaySB_omega%.1f_m%.1f_beta%.1f_lam%.1f_wD%.1f_dt%.0e_L%d_M%d_%s.out', ...
    omega, mass, beta, lambda_D, omega_D, dt, L_max, M_max, name);
    outfile = fullfile(outfolder, filename);
    fileID = fopen(outfile, 'w');

    % write down the parameters used in the simulation
    fprintf(fileID, '## Spin Boson Model Simulation Results using Tom Fays Code\n');
    fprintf(fileID, '### Parameters:\n');
    fprintf(fileID, '# omega = %f\n', omega);
    fprintf(fileID, '# m = %f\n', mass);
    fprintf(fileID, '# beta = %f\n', beta);
    fprintf(fileID, '# lambda_D = %f\n', lambda_D);
    fprintf(fileID, '# omega_D = %f\n', omega_D);
    fprintf(fileID, '# dt = %f\n', dt);
    fprintf(fileID, '# n_steps = %d\n', n_steps);
    fprintf(fileID, '# krylov_dim = %d\n', krylov_dim);
    fprintf(fileID, '# krylov_tol = %e\n', krylov_tol);
    fprintf(fileID, '# L_max = %d\n', L_max);
    fprintf(fileID, '# M_max = %d\n', M_max);
    fprintf(fileID, '# \n');
    fprintf(fileID, '## Results:\n');
    % Write the header for the results
    header = {'# t','\left<\sigma_x(t)\right>','\left<\sigma_y(t)\right>','\left<\sigma_z(t)\right>'};
    % Write the header to the file
    fprintf(fileID, '%s\t', header{1:end-1});  % Print all except last label
    fprintf(fileID, '%s\n', header{end});

    % Loop through each row and write it to the file
    for j = 1:size(O_t, 2)
        fprintf(fileID, '%f\t', t(j));
    for fxn = 1:size(O_t, 1)
        fprintf(fileID, '%f\t', O_t(fxn, j));
    end
    fprintf(fileID, '\n');  % Move to the next line after each row
    end
    fclose(fileID);
end