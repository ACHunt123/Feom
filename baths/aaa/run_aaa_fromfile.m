function run_aaa_fromfile(K, location,filename,extension,terminate,fit_wRg)
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
    if nargin < 6
        fit_wRg = false; % if fit_wRg is true, we are in the cryogenic regime (low temperature)
    end
    doplot=false;

    % Load data: x Re(f) Im(f)
    data = readmatrix(filename, 'FileType', 'text', 'Delimiter', ' ');
    x = data(:,1);
    freal = data(:,2);
    fimag = data(:,3);
    F = freal + 1i * fimag;
    Z = x;

    %%% If temperature really low, multiply by x to get rid of the singularity at x=0
    if fit_wRg 
        F=F.*Z;
    end

    %%% if terminate is true, we are allowing as many poles as possible
    if terminate
        tol=1e-10;
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
    % calculate k
    r_from_res = zeros(size(Z));  % Initialize the result vector
    for j = 1:length(res)
        r_from_res = r_from_res + res(j) ./ (Z - pol(j));
    end
    konstant = mean(F-r_from_res); %%% calculate the constant term to add back in
    r_from_res = r_from_res + konstant;  % Add the constant term
    F_noKonst = F - konstant;  % The function without the constant term
    % get half of the poles and residues (in the upper half plane)
    mask2shift = imag(pol) > 1e-5; %%% mask to select only the poles with positive imaginary part
    pol_pos = pol(mask2shift);  % Select only the poles with positive imaginary part
    res_pos = res(mask2shift);  % Select only the residues with positive imaginary part
    %%% project the error from changing to imaginary-only poles 
    if fit_wRg
        w_i=imag(pol_pos);          % pure imaginary poles
        gam_i =2*real(res_pos);     % pure real residues
        r_im_pols_real_res = zeros(size(Z));  % calculate the new function with imag poles, real residues
        for j = 1:length(res_pos)
            r_im_pols_real_res = r_im_pols_real_res + gam_i(j)*Z ./ (Z.^2 + w_i(j)^2);
        end
        diff = F_noKonst - r_im_pols_real_res;  % The error from using imag poles, real residues
        error_im_pols = sum(abs(diff).^2);
        fprintf('Error from using pure imaginary poles + real residues : %.2e\n',error_im_pols);
        %%% calculate the basis functions
        phi = zeros(length(Z), length(pol_pos));
        for j = 1:length(pol_pos)
            phi(:,j) = Z./ (Z.^2 +  w_i(j)^2);
        end
        r_im_pols = r_im_pols_real_res;
        konstant=0; %%% the constant term must be zero by symmety
    else
        
        gam_i = -2*imag(pol_pos).*imag(res_pos); %%% the new residues for imaginary-only poles
        w_i = imag(pol_pos); %%% the new poles for imaginary-only poles
        r_im_pol_im_res = zeros(size(Z)); % calculate the new function with imag poles, imag residues
        for j = 1:length(res_pos)
            r_im_pol_im_res = r_im_pol_im_res + gam_i(j) ./ (Z.^2 + w_i(j)^2);
        end
        diff = F_noKonst - r_im_pol_im_res;  % The error from using imaginary poles + imaginary residues
        error_im_pols = sum(abs(diff).^2);
        fprintf('Error from using pure imaginary poles + residues : %.2e\n', error_im_pols);
        %%% calculate the basis functions 
        phi = zeros(length(Z), length(pol_pos));
        for j = 1:length(pol_pos)
            phi(:,j) = 1./ (Z.^2 +  w_i(j)^2);
        end
        r_im_pols = r_im_pol_im_res;
    end
    %%% project the error onto the basis functions
    coeffs = phi\diff;  % Least-squares projection
    %%% add the correction to the residues
    gam_i=gam_i+coeffs;
    % make sure gam_i and w_i are real
    gam_i=real(gam_i);
    w_i=real(w_i);
    %%% calculate the new rational approximant with imaginary-only poles and the correction
    r_im_poles_corrected = zeros(size(Z)); % the corrected result vector (independent of either fit_wRg or not as the basis functions will be different)
    for j = 1:length(pol_pos)
        r_im_poles_corrected = r_im_poles_corrected + gam_i(j).* phi(:,j);
    end
    r_im_poles_corrected = r_im_poles_corrected + konstant;  % Add the constant term (0 for fit_wRg, non-zero otherwise)
    diff = F - r_im_poles_corrected;  % The error from using imaginary-only poles
    error_im_pols_corrected = sum(abs(diff).^2);
    fprintf('Error after projective correction : %.2e\n', error_im_pols_corrected);

    diff = F - r(Z);  % The error from the original approximant
    error_fullAAA = sum(abs(diff).^2);
    fprintf('Error from original AAA approximant : %.2e\n',error_fullAAA);



    if doplot
        if fit_wRg
            fac = 1; %%% factor to multiply by to get back to the original function
        else
            fac = Z; %%% factor to multiply by to get back to the original function
        end
        yy = r(Z); 
        %%% Plot the results
        plot(Z, F.*fac, 'k-', 'LineWidth', 1.5,'DisplayName',sprintf('Exact'));hold on; 
        F=0*F; %%% shift to zero mean
        plot(Z, yy.*fac-F, 'r--', 'LineWidth', 1.5, 'DisplayName',sprintf('AAA Rational Approximant (error = %.2e)', error_fullAAA)); %hold on
        % plot(Z, r_from_res-F, 'r--', 'LineWidth', 1.5, 'DisplayName',sprintf('AAA Rational Approximant from residues ')); 
        plot(Z, r_im_poles_corrected.*fac-F, 'g--', 'LineWidth', 1.5,'DisplayName',sprintf('imag-only poles + projective correction (error = %.2e)', error_im_pols_corrected));
        % Lock the current y-axis limits
        ylim_current = ylim;
        ylim manual
        xlim([-100 100])
        % Plot the red dashed line (but it won't change the y-limits)
        plot(Z, r_im_pols.*fac-F, 'b--', 'LineWidth', 1.5,'DisplayName',sprintf('imag-only poles (error = %.2e)', error_im_pols));
        % Restore original y-limits (optional if you want them unchanged)
        ylim(ylim_current);

        legend('Location', 'best');


        % legend('Original', 'AAA Rational Approximant', 'AAA Rational Approximant with imag-only poles', 'AAA Rational Approximant with imag-only poles 2')
        xlabel('x'); ylabel('f(x)')
        % title('errors in AAA Approximation of f(x)')
        title('AAA Approximation or R(w), multiplied by w')
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
    writematrix(real(konstant), fullfile(location, ['k', extension]), 'Delimiter', 'tab');
    % corrected (imag-only poles + projective correction)
    writematrix(w_i, fullfile(location, ['w_i', extension]), 'Delimiter', 'tab');
    writematrix(gam_i, fullfile(location, ['gam_i', extension]), 'Delimiter', 'tab');

    fprintf('Files saved successfully.\n');
end
