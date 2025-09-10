%  Runs the AAA algorithm on a specific function and saves results
%  to text files for further analysis.

% the function is defined in the file Sbeta.m


Z = linspace(-1000, 1000, 10000)';        % Sample points (column vector)
Z = linspace(-1000, 1000, 20000)';        % Sample points (column vector)
F = Sbeta(Z);                       % Function values at those points

tol = 3e-8;
[r, pol, res, zer, z, f, w, errvec] = aaa_algo(F, Z,tol);

xx = linspace(-1000, 1000, 20000);       % Evaluation points
yy = r(xx);                         % Evaluate the rational approximant

plot(xx, Sbeta(xx), 'k-', 'LineWidth', 1.5); hold on
plot(xx, yy, 'r--', 'LineWidth', 1.5)
legend('Original', 'AAA Rational Approximant')
xlabel('x'); ylabel('f(x)')
title('AAA Approximation of f(x)')
grid on


% % Set your threshold (e.g., 1e-5)
% residue_threshold = 1e-3;

% % Find indices where residues are "significant"
% mask = abs(res) > residue_threshold;

% % Apply mask
% pol = pol(mask);
% res = res(mask);



% Calculate the function directly from the residues
r_from_res = zeros(size(xx));  % Initialize the result vector
for j = 1:length(res)
    r_from_res = r_from_res + res(j) ./ (xx - pol(j));
end
% plot(xx, real(r_from_res), 'b-.', 'LineWidth', 1.5)      % From residues

% Calculate for just the imaginary parts
r_from_res = zeros(size(xx));  % Initialize the result vector
for j = 1:length(res)
    r_from_res = r_from_res + 1i*imag(res(j)) ./ (xx - 1i*imag(pol(j)));
    % r_from_res = r_from_res + real(res(j)) ./ (xx - 1i*imag(pol(j)));
end
plot(xx, real(r_from_res), 'b-.', 'LineWidth', 1.5)      % From residues with just the imaginary part

% r_from_mats = zeros(size(xx));  % Initialize the result vector


% legend('Original', 'AAA Rational Approximant', 'From Residues','Just imaginary part residues')
legend('Original', 'AAA Rational Approximant','Just imaginary part residues')

% Save all results as text files with tab delimiter for easy np.loadtxt reading
writematrix(real(pol), '.files/pol_real.txt', 'Delimiter', 'tab');
writematrix(imag(pol), '.files/pol_imag.txt', 'Delimiter', 'tab');

writematrix(real(res), '.files/res_real.txt', 'Delimiter', 'tab');
writematrix(imag(res), '.files/res_imag.txt', 'Delimiter', 'tab');

writematrix(real(zer), '.files/zer_real.txt', 'Delimiter', 'tab');
writematrix(imag(zer), '.files/zer_imag.txt', 'Delimiter', 'tab');

writematrix(errvec, '.files/errvec.txt', 'Delimiter', 'tab');
