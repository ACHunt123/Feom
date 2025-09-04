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
    doplot=true;

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

        % %%% Find a an upper-bound for tolerance that will give K and K+1 poles
        % tol = max_tol;          % Initial AAA tolerance
        % for KK = [K, K+1] 
        %     stepfactor = 25;          % Reset step factor for each KK
        %     % we do not reset the tolerance after thefirst iteration as it will be lower for K+1 than for K

        %     fprintf('Finding upper bound for tolerance for K=%d poles...\n', KK);
        %     while true
        %         [r, pol, res, zer, ~, ~, ~, errvec] = aaa_algo(F, Z, tol); %%% Run the AAA algorithm, cleaning up the poles afterwards
        %         pol_clean = pol(imag(pol) > 1e-10);
        %         fprintf('tol = %.1e -> %d significant poles\n', tol, numel(pol_clean));

        %         if numel(pol_clean) == KK        %%% Number of poles is desired
        %             if KK == K+1
        %                 min_tol = tol;          % Set min tolerance to current value [as this is not necessarily the best tolerance for the given K+1]
        %             elseif KK == K
        %                 max_tol = tol;          % Set max tolerance to current value [as this is not necessarily the best tolerance for the given K]
        %             end
        %             break;  

        %         elseif numel(pol_clean) > KK     %%% We have overshot
        %             tol = tol*stepfactor;          % increase tolerance back to previous value
        %             stepfactor = max(stepfactor / 10, 1.5);  % Reduce step factor BUT NOT THAT IT MUST ALWAYS > 1

        %         elseif tol < min_tol  % Too many poles, stop if tolerance is too high
        %             fprintf('Warning: Too many poles found with tolerance %.1e. Stopping.\n', tol);
        %             return;
        %         elseif tol > 100  % Too many poles, stop if tolerance is too high
        %             fprintf('Warning: too high');
        %             return;
        %         else  
        %             tol = tol / stepfactor;         % decrease tolerance to find more poles
        %         end
        %     end
        % end

        % fprintf('Found max tolerances: K=%d -> %.1e, K+1=%d -> %.1e\n', K, max_tol,K+1, min_tol);

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
                    fprintf('Press any key to continue or Ctrl+C to stop.\n');
                    pause;
                    return;
                end
                fprintf('Converged to tolerance %.5e with window %.5e\n', max_tol, tol_err);
                break;
            end            
        end
    end
    fprintf('AAA decomposition complete.\n');
    if doplot
        %%% Calculate the rational approximant with the final tolerance
        xx=Z;
        yy = r(xx);  
        %%% Calculate the rational approximant with the final tolerance, and with imaginary-only poles
        r_from_res = zeros(size(xx));
        for j = 1:length(res)
            r_from_res = r_from_res + 1i*imag(res(j)) ./ (xx - 1i*imag(pol(j)));
        end
        %%% Plot the results
        plot(x, freal, 'k-', 'LineWidth', 1.5); hold on
        plot(xx, yy, 'r--', 'LineWidth', 1.5); hold on
        plot(xx, real(r_from_res)+r(1000000), 'b-', 'LineWidth', 1.5);
        legend('Original', 'AAA Rational Approximant', 'AAA Rational Approximant with imag-only poles')
        xlabel('x'); ylabel('f(x)')
        title('AAA Approximation of f(x)')
        grid on   
        waitfor(gcf);  % Wait for user to close plot
    end
    % Save output
    k = r(1000000);  % Evaluate approximant far out for the constant term
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
