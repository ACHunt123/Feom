%%% Run AAA from file on which the function Sbeta is defined
function run_aaa_fromfile(nmodes,tol)

if nargin < 1
    nmodes = 3;  % Default number of modes for AAA expansion
end
if nargin < 2
    tol = 1e-6;  % Default tolerance for AAA algorithm
end


filename = '.files/aaa_data.txt';

% Read the data file (assuming TSV: x    f(x))
data = readmatrix(filename, 'FileType', 'text', 'Delimiter', ' ');
x = data(:,1);
freal = data(:,2);
fimag = data(:,3);

F = freal + 1i * fimag;  % Combine real and imaginary parts
Z = x;  % Use x as the sample points

[r, pol, res, zer, z, f, w, errvec] = aaa_algo(F, Z, tol, nmodes);


% Save all results as text files with tab delimiter for easy np.loadtxt reading
writematrix(real(pol), '.files/pol_real.txt', 'Delimiter', 'tab');
writematrix(imag(pol), '.files/pol_imag.txt', 'Delimiter', 'tab');

writematrix(real(res), '.files/res_real.txt', 'Delimiter', 'tab');
writematrix(imag(res), '.files/res_imag.txt', 'Delimiter', 'tab');

writematrix(real(zer), '.files/zer_real.txt', 'Delimiter', 'tab');
writematrix(imag(zer), '.files/zer_imag.txt', 'Delimiter', 'tab');

writematrix(errvec, '.files/errvec.txt', 'Delimiter', 'tab');
end