function [P,Q] = getshifts_fismith(L,A,B,U,S,V,varargin)
%
% getshifts_fismith(L,A,B,U,S,V,tol)
% 
% computes a collection of near-optimal ADI shift parameters for solving 
% AX - XB = U*S*V' using fiadi. The fiadi solution
% computed with P and Q satisfies ||X_approx||_2 \approx ||X||*tol.
% If tol is not specified, it is set to machine precision. 
%
% getshifts_fiadi(I,A,B,U,S,V, N)
% If N is vector of values, then N(i,:) optimal shift parameters for each
% AX_i - X_iB = F_i:
%
% In fiadi, the equation AX-XB = F = U*S*V' is split into the following N equations: 
%
%   AX_i - X_iB = F_i, with X_1 + ...+X_N = X, and F_i = U(:,i)*S(i,i)*V(:,i)'. 
%
% The nonzero entries of P(i,:), Q(i,:) are the optimal ADI shift parameters 
% for solving the above eqn. 
%
% See also: getshifts_smith, getshifts_fiadi.
%
% See Examples_lowrankadi34.m for more help. 
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
% Jan. 2018

%%
% Parse input: 

tolmode = 0; 
if isempty(varargin)
    tol = eps; 
    tolmode = 1; 
elseif isscalar(varargin{1})
    tol = varargin{1}; 
    tolmode = 1; 
else            % vector of # adi steps per rhs is provided
    tolmode = 0; 
    N = varargin{1};
end

% Compute Smith/adi shift parameter
% shift and rotate to real axis
s = diag(S);
r = length(s); 
 
o1 = L(1); r1 = L(2); o2 = L(3); r2 = L(4); 
T_rot = @(z) exp(-1i*angle(o1-o2))*(z-(o1+o2)/2);
T_rot_inv = @(z) z*(exp(1i*angle(o1-o2)))+ (o1+o2)/2;
o1_rot = real(T_rot(o1)); 
o2_rot = real(T_rot(o2));

% Mobius transform to symmetric disks
a = o1_rot-r1; b = o1_rot+r1; 
c = o2_rot-r2; d = o2_rot+r2;

%check for overlap
I1 = [min(a,b) max(a,b)]; 
I2 = [min(c,d) max(c,d)];
if (( I1(1) < I2(2) && I1(2) > I2(1)) || ( I2(1) < I1(2) && I2(2) > I1(1)) )
    error('ADI:getshifts_fismith: The disks containing the eigenvalues of A and B in AX-XB = F must be disjoint.')
end

cross = abs ((c - a ) *( d - b) /(c - b) /(d - a)) ; % |cross-ratio|
gam = -1 + 2*cross + 2*sqrt(cross ^2-cross) ;

B1 = -(gam +1) *(d - c) /((gam -1) + 2*(d - c)/(b - c)) ;
A1 = -2* B1 /(b - c) - 1; D1 = B1 ;
T_mob = @ (z) ( A1 * z + (B1 - c*A1)) ./(z +(D1 - c)) ;
T_mob_inv = @(z) (-z*(D1-c) + (B1 -c*A1))./(z - A1); 

%compute shifts:
p= T_rot_inv(T_mob_inv(-sqrt(gam)));  
q = T_rot_inv(T_mob_inv(sqrt(gam))); 

%make real-valued if complex part is zero
if abs(imag(p)) < 1e-15
    p = real(p);
end
if abs(imag(q)) < 1e-15
    q = real(q);
end

% Part 2: compute the req'd # of shift parameters: 
if tolmode ==1
z0 = abs((T_mob(a)+T_mob(b))/2);
mu = (z0+sqrt(gam))./(z0-sqrt(gam));
    
%%
% # iterations to apply per block:  
Nmax = ceil( log(1/(tol))/log(mu)); %max # iterations

% compute # iterations per block:
dist = min(abs(a-d), abs(b-c)); 
nrmx = estimate_normX(L,A,B,U*sqrt(S),V*sqrt(S));
Cst = (r/(tol*dist*nrmx)); 
N = min(Nmax, ceil(log(s*Cst)/log(mu)));

N(N<0) =0 ; %if the Cst*s is below tol,  do 0 ADI steps.
end

Nblk = length(N); 
P = zeros(Nblk, Nmax); % total # blocks by max # iterations.                                          
Q = P;

for j = 1:Nblk
    Ns = N(j);
    P(j,1:Ns) = p*ones(Ns,1);
    Q(j,1:Ns) = q*ones(Ns,1);
end

end
function nmx = estimate_normX(I,A,B,U,V) 
%This function applies a few iterations of fADI to estimate ||X||_2, 
% where AX - XB = U*V'.

%first estimate RHS with a rank 1 function:
U = U(:,1); 
V = V(:,1); 

%now apply a few iterations of fADI
[p,q] = getshifts_smith(I, 2); 
[Z, D, Y] = fadi(A, B, U, V, p, q); 

%estimate norm:
nmx = norm(Z)*max(abs(diag(D)))*norm(Y); 
end
