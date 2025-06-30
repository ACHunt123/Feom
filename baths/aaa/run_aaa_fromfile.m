%%% Run AAA from file on which the function Sbeta is defined

% Get the input filename from the base workspace
% filename = evalin('base', 'filename');
% % FreePoles_Jw_fromFig1.txt
filename = 'FreePoles_Jw_fromFig1.txt';

% Read the data file (assuming TSV: x    f(x))
data = readmatrix(filename, 'FileType', 'text', 'Delimiter', '\t');
x = data(:,1);
freal = data(:,2);
fimag = data(:,3);

F = freal + 1i * fimag;  % Combine real and imaginary parts
Z = x;  % Use x as the sample points

[r, pol, res, zer, z, f, w, errvec] = aaa_algo(F, Z);

% xx = linspace(-1000, 1000, 2000);       % Evaluation points
xx = Z;       % Evaluation points
yy = r(xx);                         % Evaluate the rational approximant

plot(Z, F, 'k-', 'LineWidth', 1.5); hold on
plot(xx, yy, 'r--', 'LineWidth', 1.5)
legend('Original', 'AAA Rational Approximant')
xlabel('x'); ylabel('f(x)')
title('AAA Approximation of f(x)')
grid on

% Calculate the function directly from the residues
r_from_res = zeros(size(xx));  % Initialize the result vector
for j = 1:length(res)
    r_from_res = r_from_res + res(j) ./ (xx - pol(j));
end
plot(xx, real(r_from_res), 'b-.', 'LineWidth', 1.5)      % From residues

% Save all results as text files with tab delimiter for easy np.loadtxt reading
writematrix(real(pol), '.files/pol_real.txt', 'Delimiter', 'tab');
writematrix(imag(pol), '.files/pol_imag.txt', 'Delimiter', 'tab');

writematrix(real(res), '.files/res_real.txt', 'Delimiter', 'tab');
writematrix(imag(res), '.files/res_imag.txt', 'Delimiter', 'tab');

writematrix(real(zer), '.files/zer_real.txt', 'Delimiter', 'tab');
writematrix(imag(zer), '.files/zer_imag.txt', 'Delimiter', 'tab');

writematrix(errvec, '.files/errvec.txt', 'Delimiter', 'tab');
