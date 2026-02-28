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
% system hamiltonian parameters
epsilon = 1.0 ;
Delta = 2.0 ;
% bath parameters
% beta = 1.5 ; 
beta = 10 ;
% UBO bath parameters Omega > gamma/2
Omega_UBO = 1.0 ;
gamma_UBO = 0.2 ;
lambda_UBO = 0.5 ; 

% dynamics information - parameters for the Short-Iterative Arnoldi
% Integrator
dt = 1e-2 ;
n_steps = 1000 ;
krylov_dim = 8 ;
krylov_tol = 1e-8 ;
% parmeters for heirarchy truncation using L/M truncation
% M_max = 5 ;
M_max = 3; 


% matrices of system observable operators to be returned, 
% sigma_z, sigma_x, sigma_y, and 1
O_sys = {[[1,0];[0,-1]],[[0,1];[1,0]],[[0,-1.0i];[1.0i,0]],eye(2)} ;

% initial state of the system
rho_0_sys = [[1,0];[0,0]] ;

% two objects are supplied to the HEOM dynamics function:
% "full_system" specifies the full Hamiltonian (system + bath) and the
% temperature.
% "heom_dynamics" specifies the HEOM truncation, integrator for the
% dynamics and the total propagation time and observables to be calculated.

% the full_system object contains all information about the Hamiltonian of
% the full open quantum system
full_system = struct ;
% H_sys contains the system Hamiltonian
full_system.H_sys = [[epsilon,Delta];
                     [Delta,-epsilon]];

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
% L_max_list={4,5,3,2};
% full_system_bathslist={struct("V",[[1,0];[0,-1]],"spectral_density","debye","omega_D",omega_D,"lambda_D",lambda_D),        ...  
% struct("V",[[1,0];[0,-1]],"spectral_density","debye","omega_D",omega_D,"lambda_D",lambda_D),        ...                       
% struct("V",[[1,0];[0,-1]],"spectral_density","debye (pade)","omega_D",omega_D,"lambda_D",lambda_D,"approximant_type","[N-1/N]","N_pade",M_max), ...
% struct("V",[[1,0];[0,-1]],"spectral_density","debye (pade)","omega_D",omega_D,"lambda_D",lambda_D,"approximant_type","[N-1/N]","N_pade",M_max)};
% termination_list={"low temp correction","low temp correction NZ2","none","low temp correction NZ2"};
% namelist={'ITlowtemp','NZ2lowtemp','Pade[N-1oN]','Pade[N-1_N]NZ2'};
L_max_list={4};
full_system_bathslist={struct("V",[[1,0];[0,-1]],"spectral_density","UBO","Omega",Omega_UBO,"lambda",lambda_UBO,"gamma",gamma_UBO)};
% termination_list={"low temp correction"};
termination_list={"none"};
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
    filename = sprintf('UBO_eps%.1f_D%.1f_beta%.1f_lam%.1f_Om%.1f_gam%.1f_dt%.0e_L%d_M%d_%s.out', ...
    epsilon, Delta, beta, lambda_UBO, Omega_UBO, gamma_UBO, dt, L_max, M_max, name);
    outfile = fullfile(outfolder, filename);
    fileID = fopen(outfile, 'w');

    % write down the parameters used in the simulation
    fprintf(fileID, '## Underdamped Brownian Oscillator (UBO) Simulation Results\n');
    fprintf(fileID, '### Parameters:\n');
    fprintf(fileID, '# epsilon = %f\n', epsilon);
    fprintf(fileID, '# Delta = %f\n', Delta);
    fprintf(fileID, '# beta = %f\n', beta);
    fprintf(fileID, '# lambda_UBO = %f\n', lambda_UBO);
    fprintf(fileID, '# Omega_UBO = %f\n', Omega_UBO);
    fprintf(fileID, '# gamma_UBO = %f\n', gamma_UBO); % Added gamma
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