function [H0, pos_mat, x_arr, psi_ns, E_ns] = ho_matrices(m, omega, hbar, xmin, dx, nx, ns)
    %----------------------------------------------
    % Harmonic oscillator eigenstates and matrices
    %----------------------------------------------
    % Inputs:
    %   m, omega, hbar  - mass, frequency, Planck constant
    %   xmin, dx, nx    - grid definition
    %   ns              - number of states
    %
    % Outputs:
    %   H0              - Bare hamiltonian in energy basis
    %   pos_mat         - position matrix in energy basis
    %   x_arr, psi_ns, E_ns
    %----------------------------------------------

    % Grid
    x_arr = dx * ((1:nx) + floor(xmin/dx));

    % Preallocate
    psi_ns = zeros(nx, ns);
    E_ns = zeros(1, ns);

    % Eigenstates (Hermite functions)
    alpha = m * omega / hbar;
    for i = 0:ns-1
        prefac = 1 / sqrt(2^i * factorial(i)) * (alpha/pi)^0.25;
        Hi = hermiteH(i, sqrt(alpha)*x_arr);
        psi_ns(:, i+1) = prefac .* Hi .* exp(-alpha*x_arr.^2/2);
        E_ns(i+1) = hbar * omega * (i + 0.5);
    end

    % Normalize eigenstates (optional check)
    % for i = 1:ns
    %     for j = 1:ns
    %         fprintf('Overlap(%d,%d) = %.3f\n', i, j, sum(conj(psi_ns(:,i)).*psi_ns(:,j))*dx);
    %     end
    % end

    % Bare Hamiltonian (diagonal)
    H0 = diag(E_ns);

    % Position operator in energy basis
    pos_mat = zeros(ns, ns);
    for m_ = 1:ns
        for n_ = 1:ns
            pos_mat(m_, n_) = sum(x_arr(:) .* conj(psi_ns(:,m_)) .* psi_ns(:,n_)) * dx;
        end
    end


end
