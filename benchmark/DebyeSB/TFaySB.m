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
beta = 1.0 ; 
% debye bath parameters
lambda_D = 1 ;
omega_D = 1.0 ;

% dynamics information - parameters for the Short-Iterative Arnoldi
% Integrator
dt = 1e-2 ;
n_steps = 1000 ;
krylov_dim = 8 ;
krylov_tol = 1e-8 ;
% parmeters for heirarchy truncation using L/M truncation
L_max = 5 ; 
M_max = 5 ;


% matrices of system observable operators to be returned, sigma_x, sigma_y
% sigma_z, and 1
O_sys = {[[0,1];[1,0]],[[0,-1.0i];[1.0i,0]],[[1,0];[0,-1]],eye(2)} ;

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
% baths is a cell array of structs describign each bath
full_system.baths = {struct("V",[[1,0];[0,-1]],...
    "spectral_density","debye","omega_D",omega_D,"lambda_D",lambda_D)} ;
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
heom_dynamics.heom_truncation.L_max = L_max ;
% heom_dynamics.heom_truncation.truncation_method = "frequency cut-off" ;
% heom_dynamics.heom_truncation.Gamma_cut = Gamma_cut ;
heom_dynamics.heom_truncation.heom_termination = "low temp correction" ;

% what system observables should be returned
heom_dynamics.observables = struct ;
heom_dynamics.observables.system = O_sys ;

% set the initial condition
heom_dynamics.rho_0_sys = rho_0_sys ;

% run the dynamics
[O_t,t] = runHEOMDynamics(full_system,heom_dynamics) ;


%%% save results to a text file
%make the filename 
% Create a descriptive filename
filename = sprintf('TFaySB_eps%.1f_D%.1f_beta%.1f_lam%.1f_wD%.1f_dt%.0e_L%d_M%d.txt', ...
    epsilon, Delta, beta, lambda_D, omega_D, dt, L_max, M_max);
outfile = fullfile(outfolder, filename);
fileID = fopen(outfile, 'w');
% write down the parameters used in the simulation
fprintf(fileID, '## Spin Boson Model Simulation Results using Tom Fays Code\n');
fprintf(fileID, '### Parameters:\n');
fprintf(fileID, '# epsilon = %f\n', epsilon);
fprintf(fileID, '# Delta = %f\n', Delta);
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
for i = 1:size(O_t, 2)
    fprintf(fileID, '%f\t', t(i));
for fxn = 1:size(O_t, 1)
    fprintf(fileID, '%f\t', O_t(fxn, i));
end
fprintf(fileID, '\n');  % Move to the next line after each row
end
fclose(fileID);