function [p, q] = getshifts_smith(L, varargin)
%
% computes the optimal ADI shift parameter (p,q) for solving AX - XB = U*V',
% where L = [o1 r1 o2 r2], and spectrum(A) \cup D1 = {z, |z-o1| <= r1}, 
% spectrum(B) \cup D1 = {z, |z-o2| <= r2}. 
%
% getshifts_smith(L) returns the optimal shift parameter: (p,q). 
%
% getshifts_smith(L,N) returns a vector of shift parameters: 
% P = p*ones(1,N), Q = q*ones(1,N) that can immediately be used with fadi. 
% 
% getshifts_smith(I, 'tol', Tol) determines N so that when the returned 
% vector P = p*ones(1,N), Q = q*ones(1,N) is used in fadi, the solution 
% satisfies ||X_approx||_2 < ||X||*Tol.
%
% See also: getshifts_fadi, getshifts_fismith
%
% References: 
%
% [1] Starke, Gerhard. "Near-circularity for the rational Zolotarev 
% problem in the complex plane." J. of approx. theory 70.1 (1992): 115-130.
%
% [2] Townsend, Alex, and Heather Wilber. "On the singular values of matrices 
% with high displacement rank." arXiv preprint arXiv:1712.05864 (2017).
%

% written by Heather Wilber (heatherw3521@gmail.com)
% Jan, 2018

%%
% part 1: compute ADI shift parameter: A mobius transform is used to 
% map two disks to two disks that are each symmetric about the real axis. 

% shift and rotate to real axis
o1 = L(1); r1 = L(2); o2 = L(3); r2 = L(4); 
T_rot = @(z) exp(-1i*angle(o1-o2))*(z-(o1+o2)/2);
T_rot_inv = @(z) z*(exp(1i*angle(o1-o2)))+ (o1+o2)/2;
o1_rot = real(T_rot(o1)); 
o2_rot = real(T_rot(o2));

%part 2 : Mobius transform to symmetric disks
a = o1_rot-r1; b = o1_rot+r1; 
c = o2_rot-r2; d = o2_rot+r2;  
% check for overlapping disks: 
I1 = [min(a,b) max(a,b)]; 
I2 = [min(c,d) max(c,d)]; 
if (( I1(1) < I2(2) && I1(2) > I2(1)) || ( I2(1) < I1(2) && I2(2) > I1(1)) )
    error('ADI:getshifts_smith: The disks containing the eigenvalues of A and B in AX-XB = F must be disjoint.')
end

cross = abs ((c - a) *(d - b) /(c - b)/(d - a)) ; % |cross-ratio|
gam = -1 + 2*cross + 2* sqrt (cross^2 - cross) ;

B1 = -(gam +1) *(d - c) /((gam -1) + 2*(d - c) /(b - c)) ;
A1 = -2* B1 /( b - c ) - 1;  D1 = B1 ;
T_mob = @ (z) (A1 * z + ( B1 - c * A1)) ./(z +( D1 - c)) ;
T_mob_inv = @(z) (-z*(D1-c) + (B1 - c * A1))./(z - A1); 

%compute shift parameter with inverse transform: 
p= T_rot_inv(T_mob_inv(-sqrt(gam)));   
q = T_rot_inv(T_mob_inv(sqrt(gam))); 

%make real-valued if complex part is zero
if abs(imag(p)) < 1e-15
    p = real(p);
end
if abs(imag(q)) < 1e-15
    q = real(q);
end

%%
% part 2: how many fadi iterations are needed? 

if isempty(varargin)     %just return shift parameter
    return 
elseif numel(varargin)==1 %N is specified
    N = varargin{1};  
else                     %tolerance specified
    tol =varargin{2};                   
    z0 = abs( (T_mob(a)+T_mob(b))/2 );
    mu = (z0+sqrt(gam))./(z0-sqrt(gam));    
    N = ceil( log(1/(tol))/log(mu)); 
end
    p = p*ones(1,N);
    q = q*ones(1,N);
end


