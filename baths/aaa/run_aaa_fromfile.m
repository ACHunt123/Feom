function run_aaa_fromfile(K, location,filename,extension,terminate)
    if nargin < 1
        K = 2;  % Default number of poles to find
    end
    if nargin < 2
        location = pwd;  % Default output location
    end
    if nargin < 3
        filename = 'aaa_data.txt';  % Default filename for input data
    end
    if nargin < 4
        extension = '.txt';  % Default file extension for output data
    end
    if nargin < 5
        terminate = false; % if a terminator is used, we are allowing as many poles as possible
    end
    doplot=false;

    % Load data: x Re(f) Im(f)
    data = readmatrix(filename, 'FileType', 'text', 'Delimiter', ' ');
    x = data(:,1);
    freal = data(:,2);
    fimag = data(:,3);
    F = freal + 1i * fimag;
    Z = x;

    %%% if terminate is true, we are allowing as many poles as possible
    if terminate
        tol=1e-20;
        [r, pol, res, zer, ~, ~, ~, errvec] = aaa_algo(F, Z, tol); %%% Run the AAA algorithm
        
    %%% if terminate is false, we are looking for exactly K poles    
    else    
        %%% Initialize parameters
        max_tol = 1e0;          % Stop if tolerance exceeds this
        min_tol = 1e-31;        % Stop if tolerance is less than this
        tol_err = 1e-15;    % The error window for the final tolerance SEE BINARY SEARCH 

        fprintf('Calling aaa_algo...\n');

        %%% Now do a binary search between max_tol and min_tol to find the mininum tolerance for K poles
        while true
            tol = (max_tol + min_tol) / 2;  % Start with the midpoint
            [r, pol, res, zer, ~, ~, ~, errvec] = aaa_algo(F, Z, tol); %%% Run the AAA algorithm, cleaning up the poles afterwards
            pol_clean = pol(imag(pol) > 1e-10); 
            fprintf('tol = %.1e -> %d significant poles\n', tol, numel(pol_clean)); 
            if numel(pol_clean) <= K        %%% Number of poles is desired
                max_tol = tol;          % Set max tolerance to current value [as this is not necessarily the best tolerance for the given K]
            elseif numel(pol_clean)>K
                min_tol= tol;          % Set min tolerance to current value [as this is not necessarily the best tolerance for the given K+1]
            else
                fprintf('Warning: Too few poles found with tolerance %.1e. Stopping.\n', tol);
                return;
            end
            if abs(max_tol - min_tol) < tol_err  % If the tolerances are close enough, stop
                [r, pol, res, zer, ~, ~, ~, errvec] = aaa_algo(F, Z, max_tol); %%% Run the AAA algorithm,
                pol_clean = pol(imag(pol) > 1e-10); %%% clean up poles
                % Check if we have the desired number of poles
                if numel(pol_clean) ~= K
                    fprintf('Warning: Final tolerance %.1e does not yield %d poles, found %d poles.\n', max_tol, K, numel(pol_clean));
                    % fprintf('Press any key to continue or Ctrl+C to stop.\n');
                    % pause;
                    return;
                end
                fprintf('Converged to tolerance %.5e with window %.5e\n', max_tol, tol_err);
                break;
            end            
        end
    end
    fprintf('AAA decomposition complete.\n');
    %%% project the error from changing to imaginary-only poles 
    k = r(1000000000);  % Evaluate approximant far out for the constant term
    %%% Calculate the rational approximant with the final tolerance, and with imaginary-only poles
    r_from_res = zeros(size(Z));
    for j = 1:length(res)
        r_from_res = r_from_res + 1i*imag(res(j)) ./ (Z - 1i*imag(pol(j)));
    end

    %%% get the difference from using imaginary-only poles
    imag_diff = F - r_from_res - k;  % The error from using imaginary-only poles
    error_imag = sum(abs(imag_diff).^2)


    %%% calculate the basis functions for the imaginary-only poles and project the error onto them
    mask2shift = imag(pol) > 0; %%% mask to select only the poles with positive imaginary part
    pol_pos = pol(mask2shift);  % Select only the poles with positive imaginary part
    res_pos = res(mask2shift);  % Select only the residues with positive imaginary part
    gam_i = -2*imag(pol_pos).*imag(res_pos); %%% the new residues for imaginary-only poles
    w_i = imag(pol_pos); %%% the new poles for imaginary-only poles

    phi = zeros(length(Z), length(pol_pos));
    for j = 1:length(pol_pos)
        phi(:,j) = 1./ (Z.^2 +  w_i(j)^2);
    end
    %%% project the error onto the basis functions
    coeffs = phi\imag_diff;  % Least-squares projection
    %%% add the correction to the residues
    gam_i=gam_i+coeffs;
    % make sure gam_i and w_i are real
    gam_i=real(gam_i);
    w_i=real(w_i);
    %%% calculate the new rational approximant with imaginary-only poles and the correction
    r_from_res_2 = zeros(size(Z));
    for j = 1:length(pol_pos)
        r_from_res_2 = r_from_res_2 + gam_i(j)./ (Z.^2 + w_i(j)^2);
    end
     imag_diff = F - r_from_res_2 - k;  % The error from using imaginary-only poles
    error_corrected_imag = sum(abs(imag_diff).^2)

    diff = F - r(Z);  % The error from the original approximant
    error_cmplx = sum(abs(diff).^2)


    if doplot
        yy = r(Z); 
        %%% Plot the results
        % plot(Z, freal, 'k-', 'LineWidth', 1.5); hold on
        plot(Z, yy-freal, 'r--', 'LineWidth', 1.5); hold on
        plot(Z, real(r_from_res_2)+r(1000000000)-freal, 'g-', 'LineWidth', 1.5); hold on
        % Lock the current y-axis limits
        ylim_current = ylim;
        ylim manual
        % Plot the red dashed line (but it won't change the y-limits)
        plot(Z, real(r_from_res)+r(1000000000)-freal, 'b-', 'LineWidth', 1.5); hold on
        % Restore original y-limits (optional if you want them unchanged)
        ylim(ylim_current);

        legend( 'AAA Rational Approximant','imag-only poles + projective correction', 'AAA Rational Approximant with imag-only poles')
        legend( sprintf('AAA Rational Approximant (error = %.2e)', error_cmplx), ...
                sprintf('imag-only poles + projective correction (error = %.2e)', error_corrected_imag), ...
                sprintf('AAA Rational Approximant with imag-only poles (error = %.2e)', error_imag));


        % legend('Original', 'AAA Rational Approximant', 'AAA Rational Approximant with imag-only poles', 'AAA Rational Approximant with imag-only poles 2')
        xlabel('x'); ylabel('f(x)')
        title('errors in AAA Approximation of f(x)')
        grid on   
        waitfor(gcf);  % Wait for user to close plot
    end
    % Save output
    fprintf('Saving results to folder: %s\n', location);
    % non-corrected
    writematrix(real(pol), fullfile(location, ['pol_real', extension]), 'Delimiter', 'tab');
    writematrix(imag(pol), fullfile(location, ['pol_imag', extension]), 'Delimiter', 'tab');
    writematrix(real(res), fullfile(location, ['res_real', extension]), 'Delimiter', 'tab');
    writematrix(imag(res), fullfile(location, ['res_imag', extension]), 'Delimiter', 'tab');
    writematrix(real(zer), fullfile(location, ['zer_real', extension]), 'Delimiter', 'tab');
    writematrix(imag(zer), fullfile(location, ['zer_imag', extension]), 'Delimiter', 'tab');
    writematrix(errvec, fullfile(location, ['errvec', extension]), 'Delimiter', 'tab');
    writematrix(k, fullfile(location, ['k', extension]), 'Delimiter', 'tab');
    % corrected (imag-only poles + projective correction)
    writematrix(w_i, fullfile(location, ['w_i', extension]), 'Delimiter', 'tab');
    writematrix(gam_i, fullfile(location, ['gam_i', extension]), 'Delimiter', 'tab');

    fprintf('Files saved successfully.\n');
end
