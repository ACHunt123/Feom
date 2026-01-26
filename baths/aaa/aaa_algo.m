%Taken from
%Y. Nakatsukasa, O. S`ete, and L. N. Trefethen, 
% "The AAA Algorithm for Rational Approximation," 
% SIAM J. Sci. Comput. 40, A1494–A1522 (2018).
function [r,pol,res,zer,z,f,w,errvec] = aaa_algo(F,Z,tol,mmax)
% aaa rational approximation of data F on set Z
% [r,pol,res,zer,z,f,w,errvec] = aaa(F,Z,tol,mmax)
%
% Input: F = vector of data values, or a function handle
% Z = vector of sample points
% tol = relative tolerance tol, set to 1e-13 if omitted
% mmax: max type is (mmax-1,mmax-1), set to 100 if omitted
%
% Output: r = AAA approximant to F (function handle)
% pol,res,zer = vectors of poles, residues, zeros
% z,f,w = vectors of support pts, function values, weights
% errvec = vector of errors at each step
M = length(Z); % number of sample points
if nargin<3, tol = 1e-13; end % default relative tol 1e-13
if nargin<4, mmax = 100; end % default max type (99,99)
if ~isfloat(F), F = F(Z); end % convert function handle to vector
Z = Z(:); F = F(:); % work with column vectors
SF = spdiags(F,0,M,M); % left scaling matrix
J = 1:M; z = []; f = []; C = []; % initializations
errvec = []; R = mean(F);
for m = 1:mmax % main loop
[~,j] = max(abs(F-R)); % select next support point
z = [z; Z(j)]; f = [f; F(j)]; % update support points, data values
J(J==j) = []; % update index vector
C = [C 1./(Z-Z(j))]; % next column of Cauchy matrix
Sf = diag(f); % right scaling matrix
A = SF*C - C*Sf; % Loewner matrix
[~,~,V] = svd(A(J,:),0); % SVD
w = V(:,m); % weight vector = min sing vector
N = C*(w.*f); D = C*w; % numerator and denominator
R = F; R(J) = N(J)./D(J); % rational approximation
err = norm(F-R,inf);
errvec = [errvec; err]; % max error at sample points
if err <= tol*norm(F,inf), break, end % stop if converged
end
r = @(zz) feval(@rhandle,zz,z,f,w); % AAA approximant as function handle
[pol,res,zer] = prz(r,z,f,w); % poles, residues, and zeros
% [r,pol,res,zer,z,f,w] = cleanup(r,pol,res,zer,z,f,w,Z,F); % remove Frois. doublets (optional)
function [pol,res,zer] = prz(r,z,f,w) % compute poles, residues, zeros
m = length(w); B = eye(m+1); B(1,1) = 0;
E = [0 w.'; ones(m,1) diag(z)];
pol = eig(E,B); pol = pol(~isinf(pol)); % poles
dz = 1e-5*exp(2i*pi*(1:4)/4);
res = r(bsxfun(@plus,pol,dz))*dz.'/4; % residues
E = [0 (w.*f).'; ones(m,1) diag(z)];
zer = eig(E,B); zer = zer(~isinf(zer)); % zeros
function r = rhandle(zz,z,f,w) % evaluate r at zz
zv = zz(:); % vectorize zz if necessary
CC = 1./bsxfun(@minus,zv,z.'); % Cauchy matrix
r = (CC*(w.*f))./(CC*w); % AAA approx as vector
ii = find(isnan(r)); % find values NaN = Inf/Inf if any
for j = 1:length(ii)
r(ii(j)) = f(find(zv(ii(j))==z)); % force interpolation there
end
r = reshape(r,size(zz)); % AAA approx
% save('aaa_output.mat', 'z', 'f', 'w');  % Save for Python
