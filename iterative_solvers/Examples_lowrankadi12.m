% These basic examples show how fadi (Example 1) and 
% fiadi (Example 2) can be used to find a low rank approximation to
% X in AX - XB = U*V'  when the eigenvalues of A and B are real.
% (See Examples_lowrankadi34.m for complex-valued eigenvalues).
% The comments discuss when it is appropriate to use these methods. 
%
% Heather Wilber, Jan. 2018

%%
% Example 1: Using fADI
%
% We'll start by creating A, B, and a right-hand side F = U*V'.
% fadi is highly effective when rank(F) is small, 
% A and B are normal matrices, and the eigenvalues of A and B are are contained in 
% disjoint, real intervals, or disjoint disks in the complex plane. In this 
% example, we set up A and B with eigenvalues contained in disjoint intervals.
% See Examples_lowrankadi34.m for another example. 

n = 256; 
a = 0; 
b = 5; 
c = -21; 
d = -1; 
A = rand(n,n); 
[Q, ~] = qr(A); 
dg  = linspace(a,b, n); 
A = Q*diag(dg)*Q'; % A is normal, has real eigenvalues in the interval [a b]

B = rand(n,n); 
[Q, ~] = qr(B); 
dg = linspace(c,d, n); 
B = Q*diag(dg)*Q'; % B is normal, has real eigenvalues in the interval [c d]

%%
% Usually, one also requires that A and B are sparse or structured, so that
% linear solves involving them are fast. This ensures that ADI-based methods
% are efficient. For this example, we will ignore this requirement.

% The right-hand side is a rank r matrix: 
r = 1; 
U = rand(n,r); 
V = rand(n,r); 

%%
% To use fadi to find X in AX - XB = UV', we require 
% a set of ADI shift parameters. There are many possible options 
% (for example, see [1]), and any set of parameters works with this code. 
% However, when the eigenvalues of A and B are 
% contained in disjoint intervals, we can use the following command to get 
% optimal shift parameters (see [2]and [4]):  

tol = 1e-10; 
I = [a b c d]; 
[p,q] = getshifts_adi(I, 'tol', tol); 

%%
% We can also compute exactly N shift parameters:

N = 3;
[p2, q2] = getshifts_adi(I, N); 

% You'll notice that the shift parameters aren't nested; there is no reason
% that the entries in p2 should match the first three entries of p. 

%%
% Now we have everything we need to call fadi.  

[Z, D, Y] = fadi(A, B, U, V, p, q);
%%
% We compute the solution using lyap to test the accuracy.

Xt = lyap(A, -B,-U*V'); 
norm(Xt - Z*D*Y')/norm(Xt)

[~,rnk] = size(Z) %the rank of the uncompressed solution = (number of shift parameters)*r  

%%
% It can be beneficial to compress the low rank factors when r > 1:
r = 10; 
uncompressed_rank = r*length(p)

U = rand(n, r); 
V = rand(n,r); 

[Z, D, Y] = fadi(A, B, U, V, p, q,tol);

Xt = lyap(A, -B,-U*V'); 
norm(Xt - Z*D*Y')/norm(Xt)
[~,rnk] = size(Z)

%%
% Example 2: Using fiadi
%
% fiadi is a generalization of fadi that can be used when F is of moderate
% to high rank, but has singular values that decay rapidly. 
% We now set up a right hand side F where fiadi is useful. 

r = 50; 
U = rand(n,r); 
V = rand(n,r); 
[U,~,V] = lowrank_svd(U, V); 

%replace S with rapidly decaying sequence of singular values:
rho = 3;  
s = rho.^(-(1:r)); 
S = diag(s); 
% Now we have a RHS F = U*S*V'. 

%%
% The FIADI method works by splitting 
% AX -XB = USV' into N equations, 
% 
%  AX_i - X_iB = F_i, where X_1 + ... + X_N = X, and F_i = U(:,i)*S(i,i)*V(:,i)'. 
%
% Then, fadi is applied to each equation AX_i - X_i = F_i, but with ADI shift
% parameters that account for the influence of the singular values of F. 
% Each i requires its own set of ADI shift parameters (see [5]). 
% In the case where the eigenvalues of A and B are contained in disjoint 
% real intervals, near-optimal parameters can be computed with the 
% getshifts_fiadi command: 

[P, Q] = getshifts_fiadi(I,A,B, U, S, V, tol); 

%This time, matrices P and Q containing several sets of shift parameters 
% are returned. Specifically, P(i,:) contains the shift parameters 
% for solving AX_i - X_i = F_i above. 
%
% Only the nonzero entries of P and Q are used as shift parameters. Looking 
% at P, we can see that as the singular values of S get smaller, 
% fewer shift parameters are needed: 
 
surf(abs(P)), view(2), shading flat, colorbar
xlabel('shift parameters')
ylabel('RHS singular value index (i in S(i,i))')

%%
% We can use P and Q with fiadi. The solution is compressed to the 
% relative accuracy given by tol. 

[ZZ,DD,YY] = fiadi(A,B,U,S,V,P,Q,tol); 

Xt = lyap(A, -B,-U*S*V'); 
norm(Xt - ZZ*DD*YY')/norm(Xt)


%%
% A note on fiadi vs. fadi: These two methods are very similar and their
% comparative performances depend strongly on methods of implementation. 
%  
% In general, fiadi is recommended if rank(F) is very large with 
% rapidly decaying singular values, and the influence of the small singular
% values is important. 


%%
% See Examples_lowrankadi34.m, Example_poisson_fd.m and 
% Examples_adi1.m for related discussion.

%%
% References: 
%
%   [1] Benner, Peter, Ren-Cang Li, and Ninoslav Truhar. 
%   "On the ADI method for Sylvester equations." J. of Comp. and App. Math.
%   233.4 (2009): 1035-1045.
% 
%   [2] B. Beckermann and A. Townsend, 
%   "On the singular values of matrices with displacement structure."
%   SIMAX, 38 (2017): pp. 1227-1248.
%
%   [3] Fortunato, Daniel, and Alex Townsend. 
%   "Fast Poisson solvers for spectral methods." arXiv preprint arXiv:1710.11259 (2017).
%
%   [4] Lu, An, and Eugene L. Wachspress. 
%   "Solution of Lyapunov equations by alternating direction implicit iteration." 
%   Comp. & Math. with Appl., 21.9 (1991): pp. 43-58.
% 
%   [5] Townsend, Alex, and Heather Wilber. 
%   "On the singular values of matrices with high displacement rank." 
%   arXiv preprint arXiv:1712.05864 (2017).
%


