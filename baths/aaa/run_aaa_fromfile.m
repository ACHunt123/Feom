function run_aaa_fromfile(K, location,filename,extension)
    if nargin < 1
        K = 2;  % Default number of poles to find
    end
    if nargin < 2
        location = pwd;  % Default output location
    end
    doplot=false;
    tol = 1e-20;            % Initial AAA tolerance
    stepfactor = 25;        % Tolerance multiplier
    max_tol = 1e1;         % Stop if tolerance exceeds this

    % Load data: x Re(f) Im(f)
    data = readmatrix(filename, 'FileType', 'text', 'Delimiter', ' ');
    x = data(:,1);
    freal = data(:,2);
    fimag = data(:,3);
    F = freal + 1i * fimag;
    Z = x;

    fprintf('Calling aaa_algo...\n');

    % Try increasing tolerance until we get K poles with significant imaginary part
    while true
        [r, pol, res, zer, ~, ~, ~, errvec] = aaa_algo(F, Z, tol);

        % Remove near-real poles, keep only positive imaginary poles
        pol_clean = pol(imag(pol) > 1e-10);
        
        fprintf('tol = %.1e -> %d significant poles\n', tol, numel(pol_clean));

        if numel(pol_clean) == K
            break;  % Desired number of poles found
        
        elseif numel(pol_clean) < K %we have overshot
            tol = tol/stepfactor;          % Decrease tolerance back to previous value
            stepfactor = stepfactor / 10;  % Reduce step factor

        elseif tol > max_tol  % Too many poles, stop if tolerance is too high
            fprintf('Warning: Too many poles found with tolerance %.1e. Stopping.\n', tol);
            return;
        else  % Too many poles, increase tolerance
            tol = tol * stepfactor;         % Increase tolerance
        end
    end

    fprintf('AAA decomposition complete.\n');
    if doplot
        % Plotting
        xx = linspace(-1000, 1000, 2000);
        yy = r(xx);  % Rational approximant

        plot(x, freal, 'k-', 'LineWidth', 1.5); hold on
        plot(xx, yy, 'r--', 'LineWidth', 1.5)
        legend('Original', 'AAA Rational Approximant')
        xlabel('x'); ylabel('f(x)')
        title('AAA Approximation of f(x)')
        grid on

        % From residues
        r_from_res = zeros(size(xx));
        for j = 1:length(res)
            r_from_res = r_from_res + res(j) ./ (xx - pol(j));
        end
        plot(xx, real(r_from_res), 'b-.', 'LineWidth', 1.5)

        % Imaginary-only poles and residues
        r_from_res = zeros(size(xx));
        for j = 1:length(res)
            r_from_res = r_from_res + 1i*imag(res(j)) ./ (xx - 1i*imag(pol(j)));
        end
        plot(xx, real(r_from_res), 'g-.', 'LineWidth', 1.5)

        
        waitfor(gcf);  % Wait for user to close plot
    end
    % Save output
    k = r(1000000);  % Evaluate approximant far out
    fprintf('Saving results to folder: %s\n', location);
    writematrix(real(pol), fullfile(location, ['pol_real', extension]), 'Delimiter', 'tab');
    writematrix(imag(pol), fullfile(location, ['pol_imag', extension]), 'Delimiter', 'tab');
    writematrix(real(res), fullfile(location, ['res_real', extension]), 'Delimiter', 'tab');
    writematrix(imag(res), fullfile(location, ['res_imag', extension]), 'Delimiter', 'tab');
    writematrix(real(zer), fullfile(location, ['zer_real', extension]), 'Delimiter', 'tab');
    writematrix(imag(zer), fullfile(location, ['zer_imag', extension]), 'Delimiter', 'tab');
    writematrix(errvec, fullfile(location, ['errvec', extension]), 'Delimiter', 'tab');
    writematrix(k, fullfile(location, ['k', extension]), 'Delimiter', 'tab');

    fprintf('Files saved successfully.\n');
end
