%  Runs the AAA algorithm on a specific function and saves results
%  to text files for further analysis.

% the function is defined in the file Sbeta.m


Z = linspace(-1000, 1000, 10000)';        % Sample points (column vector)
Z = linspace(-100, 100, 20000)';        % Sample points (column vector)
F = Sbeta(Z);                       % Function values at those points

tol = 1e-10;
tol = 1e-11;

[r, pol, res, zer, z, f, w, errvec] = aaa_algo(F, Z,tol);

xx = linspace(-100, 100, 20000);       % Evaluation points
yy = r(xx);                         % Evaluate the rational approximant

plot(xx, Sbeta(xx), 'k-', 'LineWidth', 1.5,'DisplayName', 'original'); hold on
plot(xx, yy, 'r--', 'LineWidth', 1.5,'DisplayName', 'AAA Rational Approximant'); 


if 0 %% Remove small residues and corresponding poles
    sprintf('removing small residues and corresponding poles');
    % Set your threshold (e.g., 1e-5)
    residue_threshold = 1e-3;
    % Find indices where residues are "significant"
    mask = abs(res) > residue_threshold;
    % Apply mask
    pol = pol(mask);
    res = res(mask);
end

% Calculate the function directly from the residues
r_from_res = zeros(size(Z));  % Initialize the result vector
for j = 1:length(res)
    r_from_res = r_from_res + res(j) ./ (Z - pol(j));
end
k = r(Z(10000))-r_from_res(10000)
plot(Z, real(r_from_res)+k, 'b-.', 'LineWidth', 1.5, 'DisplayName','calculated from residues');      % From residues

if 0    % Calculate for just the imaginary parts
    r_from_res = zeros(size(Z));  % Initialize the result vector
    for j = 1:length(res)
        r_from_res = r_from_res + 1i*imag(res(j)) ./ (Z - 1i*imag(pol(j)));
    end
    plot(Z, real(r_from_res), 'b-.', 'LineWidth', 1.5,'DisplayName', 'imag res and pols');     % From residues with just the imaginary part
end

% get half of the poles and residues (in the upper half plane)
mask2shift = imag(pol) > 0; %%% mask to select only the poles with positive imaginary part
pol_pos = pol(mask2shift);  % Select only the poles with positive imaginary part
res_pos = res(mask2shift);  % Select only the residues with positive imaginary part
w_i=imag(pol_pos);
gam_i =2*real(res_pos);
r_from_res = zeros(size(Z));  % Initialize the result vector
for j = 1:length(res_pos)
    r_from_res = r_from_res + gam_i(j)*Z ./ (Z.^2 + w_i(j)^2);
end
plot(Z, real(r_from_res), 'm-.', 'LineWidth', 1.5, 'DisplayName', 'imag res real pols');     % From residues with just the imaginary part      

plot(Z, 1./(2*abs(Z)), 'b-.', 'LineWidth', 1.5, 'DisplayName', '1/2|w|');     % From residues with just the imaginary part      


diff = F - real(r_from_res);  % The error from using imaginary-only poles
%%% calculate the basis functions for the imaginary-only poles and project the error onto them
phi = zeros(length(Z), length(pol_pos));
for j = 1:length(pol_pos)
    phi(:,j) = Z./ (Z.^2 +  w_i(j)^2);
end
%%% project the error onto the basis functions
coeffs = phi\diff;  % Least-squares projection
%%% add the correction to the residues
gam_i=gam_i+coeffs;

r_from_res = zeros(size(xx));  % Initialize the result vector
for j = 1:length(res_pos)
    r_from_res = r_from_res + gam_i(j)*xx ./ (xx.^2 + w_i(j)^2);
end
sprintf('number of poles used in the imaginary-only approximation: %d', length(res_pos))
plot(xx, real(r_from_res), 'g-.', 'LineWidth', 1.5, 'DisplayName', 'imag res real pols CORRECTED');        % From residues with just the imaginary part


xlim([-100 100])
xlabel('x'); ylabel('f(x)')
title('AAA Approximation of f(x)')
grid on
% r_from_mats = zeros(size(xx));  % Initialize the result vector


% legend('Original', 'AAA Rational Approximant', 'From Residues','Just imaginary part residues')

legend show

% Save all results as text files with tab delimiter for easy np.loadtxt reading
writematrix(real(pol), '.files/pol_real.txt', 'Delimiter', 'tab');
writematrix(imag(pol), '.files/pol_imag.txt', 'Delimiter', 'tab');

writematrix(real(res), '.files/res_real.txt', 'Delimiter', 'tab');
writematrix(imag(res), '.files/res_imag.txt', 'Delimiter', 'tab');

writematrix(real(zer), '.files/zer_real.txt', 'Delimiter', 'tab');
writematrix(imag(zer), '.files/zer_imag.txt', 'Delimiter', 'tab');

writematrix(errvec, '.files/errvec.txt', 'Delimiter', 'tab');
