function [p, q] = getshifts_adi(I, varargin)
%
% computes optimal ADI shift parameters for solving AX - XB = F
% where I = [a b c d], and spectrum(A) \cup [a, b] and spectrum(B) \cup [c, d]. 
% 
% getshifts_adi(I, N) computes N ADI shift parameters. 
% 
% getshifts_adi(I, 'tol', Tol) computes as many shift parameters 
% as needed so that when they are used with adi or fadi, 
% ||X_approx||_2 < ||X||*Tol
%
% See also: getshifts_fiadi, getshifts_smith
%
% References: 
%  [1] Lu, An, and Eugene L. Wachspress. 
%  "Solution of Lyapunov equations by alternating direction implicit iteration." 
%  Comp. & Math. with Appl., 21.9 (1991): pp. 43-58.
%
%  [2] B. Beckermann and A. Townsend, 
%  "On the singular values of matrices with displacement structure."
%  SIMAX, 38 (2017): pp. 1227-1248. 

% written by Heather Wilber (heatherw3521@gmail.com)
% Jan, 2018
%%

a = I(1); 
b = I(2); 
c = I(3); 
d = I(4);

%check if intervals are overlapping: 
I1 = [min(a,b) max(a,b)]; 
I2 = [min(c,d) max(c,d)]; 
if (( I1(1) < I2(2) && I1(2) > I2(1)) || ( I2(1) < I1(2) && I2(2) > I1(1)) )
    error('ADI:getshifts_adi: The intervals containing the eigenvalues of A and B in AX-XB = F must be disjoint.')
end

% compute mobius transform
[~, Tinv, gam, cr] = mobiusT(I); 

%if tol is specified, find N: 
if numel(varargin) > 1
    tol= varargin{2};
    N = ceil(1/pi^2*log(4/tol)*log(16*cr));
else
    N = varargin{1}; 
end
    
%we estimate elliptic integrals when 1/gam is small. 
if gam > 1e7 
    K = (2*log(2)+log(gam)) + (-1+2*log(2)+log(gam))/gam^2/4;
    m1 = 1/gam^2; 
    u = (1/2:N-1/2)*K/N; 
    dn = sech(u) + .25*m1*(sinh(u).*cosh(u)+u).*tanh(u).*sech(u);        
else
kp = 1-(1/gam)^2; 
K = ellipke(kp); 
[~, ~, dn] = ellipj((1:2:(2*N-1))*K/2/N, kp);
end

%optimal shift parameters Z(T([a b]), T([c,d])) 
p1 = gam*dn;

%solve for zeros and poles of rat function on [a b] and [c d]. 
    p = Tinv(-p1); 
    q = Tinv(p1);
end


function [T, Tinv, gam, M] = mobiusT(I)
%given I = [a b c d] where [a b] [c d] are two disjoint intervals 
% on the real line, T(I) maps to the four points [-gamma, -1, 1, gamma]. 
% M is the cross-ratio. 

a = I(1); 
b = I(2); 
c = I(3); 
d = I(4); 

%parameters
M = abs((c-a)*(d-b)/((c-b)*abs(d-a))); 
gam = -1+2*M+2*sqrt(M^2-M); 
A = -gam*a*(1-gam)+gam*(c-gam*d)+c*gam-gam*d; 
B = -gam*a*(c*gam-d)-a*(c*gam-gam*d)-gam*(c*d-gam*d*c); 
C = a*(1-gam)+gam*(c-d)+c*gam-d; 
D = -gam*a*(c-d)-a*(c-gam*d)+c*d-gam*d*c; 

T = @(z) (A*z+B)./(C*z+D);
Tinv = @(z) (D*z-B)./(-C*z+A);
end
    