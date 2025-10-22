% % this file contains the function Sbeta that we will approximate
function Sb = Sbeta(omega)
    beta = 10;
    beta = 1000000000;

    hbar = 1;
    m = 1;
    u = beta * hbar * omega / 2;

    % get J(w)
    eta=1;
    gamma=1;
    Jw = eta*gamma.*omega./(omega.^2 + gamma^2);


    % Positive omega: ω^s * exp(-ω / ω_c)
    % Jw(omega > 0) =  omega(omega > 0).^s .* exp(-omega(omega > 0) / wc);

    % Negative omega: -|ω|^s * exp(-|ω| / ω_c)
    % Jw(omega < 0) = -abs(omega(omega < 0)).^s .* exp(-abs(omega(omega < 0)) / wc);

    % Sb = omega.*(1 ./ (m * beta * omega.^2)) .* ((beta * hbar * omega / 2) .* coth(beta * hbar * omega / 2) - 1);
    Sb = (1 ./ (m * beta * omega.^2)) .* ((beta * hbar * omega / 2) .* coth(beta * hbar * omega / 2) - 1);
    
    
    % Sb = Jw.*(coth(u) - 1./u);
    % Sb = Jw.*(coth(u) + 1);  
    % y=coth(omega / 2);
    % y=coth(omega);
    % y =  omega./(1 + omega.^2-omega);

end


% function wJw = Sbeta(omega)
%     s = 1/2;         % You can adjust this
%     wc = 1.0;        % Cutoff frequency (ω_c)

%     Jw = zeros(size(omega));  % Preallocate output array

%     % Positive omega: ω^s * exp(-ω / ω_c)
%     Jw(omega > 0) =  omega(omega > 0).^s .* exp(-omega(omega > 0) / wc);

%     % Negative omega: -|ω|^s * exp(-|ω| / ω_c)
%     Jw(omega < 0) = -abs(omega(omega < 0)).^s .* exp(-abs(omega(omega < 0)) / wc);

%     wJw = Jw;%.* omega;  % Multiply by omega
% end