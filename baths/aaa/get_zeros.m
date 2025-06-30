
% Read the coefficints data file (assuming TSV: x    f(x))
filename = '.coefficients/FreePoles_Sw_minus_coeffs.txt';
data = readmatrix(filename, 'FileType', 'text', 'Delimiter', '\t');
coeffs = fliplr(data(:,1));

% read the poles and residues from another file
filename = '.coefficients/FreePoles_Sw_minus_poles.txt';
data = readmatrix(filename, 'FileType', 'text', 'Delimiter', '\t');
residues = data(:,1) + 1i * data(:,2);  % Combine real and imaginary parts
poles = data(:,3) + 1i * data(:,4);  % Combine real and imaginary parts


% Calculate the roots of the polynomial defined by the coefficients
r=roots(coeffs);



% test that the roots are correct
% define function to evaluate the polynomial at the roots from the poles and residues
% and check if they are close to zero
tolerance = 1; % tolerance for checking if close to zero
valid_roots = [];
for i = 1:length(r)

    value =0;  % Initialize the result vector
    for j = 1:length(residues)
        value = value + residues(j) ./ (r(i) - poles(j));
    end

    if abs(value) > tolerance
        fprintf('Root %d is not a root of the polynomial: %e\n', i, value);
    else
        fprintf('Root %d is a valid root of the polynomial.\n', i);
        valid_roots(end+1) = r(i);  % append valid root
    end
end

% Do a plot of the roots, only if they are valid
figure;
plot(real(valid_roots), imag(valid_roots), 'o');
xlabel('Real Part');
xlim([-5 5]);
ylabel('Imaginary Part');
title('Roots of the Polynomial');
grid on;  