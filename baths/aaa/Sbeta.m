% this file contains the function Sbeta that we will approximate
function y = Sbeta(omega)
    beta = 1;
    hbar = 1;
    m = 1;
    u = beta * hbar * omega / 2;
    y = (1 ./ (m * beta * omega.^2)) .* ((beta * hbar * omega / 2) .* coth(beta * hbar * omega / 2) - 1);
    % y=coth(omega / 2);
    % y=coth(omega);
    % y =  omega./(1 + omega.^2-omega);

end