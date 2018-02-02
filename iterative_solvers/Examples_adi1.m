% This example shows how to use the adi function to solve AX-XB = F.
% 
% This code is more efficient that lyap in the following situation:
%
% (i) A and B are large and sparse so that fast linear solves with A and B 
% are possible. 
% (ii) A and B are normal matrices
% (iii) the eigenvalues of A and B are either contained in real disjoint 
% intervals, or disjoint disks in the complex plane. See [1,2,3,4].
%
% One can still use the code when (ii) and (iii) are not satisfied, but
% it will not be as efficient. (TO DO: add getshifts_penzl and demo this)
% 
% ADI can also be used to solve AX -XB = F in low rank form, see fadi and
% fiadi, and Examples_lowrank12.m, Examples_lowrank34.m for details.

% written by heather wilber (heatherw3521@gmail.com)


%%
% To use adi, we first need a relevant problem of the form AX-XB = F.
% A classic example is given in Examples_poisson_fd.m. 
% We will ignore requirement (i) in this example, since it isn't important
% for functionality. We begin by constructing A and B to have real
% eigenvalues in disjoint intervals [a b] and [c d], respectively. 

clear all 
n = 100; 
a = -100; b = -.01; c = .01; d = 20; 
A = rand(n,n);
[Q, ~] = qr(A); 
D = linspace(a,b,n); 
A = Q*diag(D)*Q'; % A is normal with eigenvalues in [a b]

B = rand(n,n); 
[Q, ~] = qr(B); 
D = linspace(c,d,n); 
B = Q*diag(D)*Q'; % B is normal with eigenvalues in [c d]. 

F = rand(n,n); 

%%
% ADI is an iterative algorithm that solves AX -XB = F. At each iteration i, 
% a pair of numbers (p_i, q_i) called a shift parameter is required. 
% The effectiveness of ADI depends strongly on these parameters. When it is
% known that the eigenvalues of A and B are contained in intervals [a b], [c d],
% respectively, then optimal shift parameters can be computed (See [1,2]).

% The getshifts_adi function returns these parameters. Here, we compute the
% first 5: 

I = [a b c d]; 
N = 5;
[p5, q5] = getshifts_adi(I, N); 

%%
% let's try using these shift parameters to solve for X with adi: 

X = adi(A, B, F, p5, q5); 

% We will compare our solution to the solution found using lyap: 

Xt = lyap(A, -B, -F); 

norm(X - Xt)/norm(Xt)

%%
% Using more than 5 shift parameters (i.e., more iterations) results in 
% better accuracy. We can specify a desired relative accuracy and compute
% exactly the number of shifts we need: 
tol = 1e-10; 
[p, q] = getshifts_adi(I, 'tol', tol); 

% Note that shift parameters are not 'nested'; the entries in p5 do not 
% match the first 5 entries in p. Let's try using these parameters in adi: 

X = adi(A, B, F, p, q); 

norm(X - Xt)/norm(Xt)

%%
% Optimal shift parameters can also be computed when the eigenvalues of 
% A and B are each contained in two disjoint disks in the complex plane
% (See [3,4]). We now set up A and B in this way: 

%%
% we choose A to have eigenvalues in a disk with center o1  and
% radius r1, and B to have eigenvalues in a disk with center o2, radius r2.

o1 = -10+0i; 
r1 = 5; 
o2 = 2 + 5i; 
r2 = 3; 

A = rand(n, n); 
[Q, ~] = qr(A,0); 
t = -pi + rand(n,1)*2*pi; 
rad = rand(n,1);
eigsA = rad*(r1).*exp(1i*t)+o1;
A = Q*diag(eigsA)*Q'; % A is normal and has the desired eigenvalues

B = rand(n, n);
[Q, ~] = qr(B, 0); 
t = -pi + rand(n,1)*2*pi; 
rad = rand(n,1);
eigsB = rad*(r2).*exp(1i*t)+o2;
B = Q*diag(eigsB)*Q'; % B is normal and has the desired eigenvalues

%%
% In this case, the shift parameter (p_i, q_i) is the same 
% at every iteration i. The adi method is therefore equivalent to Smith's 
% method. We can compute this parameter as follows: 

I = [o1 r1 o2 r2];
[p,q] = getshifts_smith(I); 

%%
% To use the parameter in adi and perform N iterations, we need
% input vectors of the form P = p*ones(1,N), Q = q*ones(1,N); 
% If we specify N in getshifts_smith, this is done automatically: 

[P, Q] = getshifts_smith(I,N); 

%%
% As before, we can also get exactly the right number of shifts to 
% achieve a desired relative accuracy: 

[P, Q] = getshifts_smith(I, 'tol', tol); 

X = adi(A, B, F, P, Q); 

Xt = lyap(A, -B, -F); 

norm(Xt - X)/norm(Xt)

%%
% See also Examples_poisson_fd.m, Examples_lowrank12.m,
% Examples_lowrank34.m
%
%% References: 
%
%   [1] B. Beckermann and A. Townsend, 
%   "On the singular values of matrices with displacement structure."
%   SIMAX, 38 (2017): pp. 1227-1248.
%
%   [2] Lu, An, and Eugene L. Wachspress. 
%   "Solution of Lyapunov equations by alternating direction implicit iteration." 
%   Comp. & Math. with Appl., 21.9 (1991): pp. 43-58.
% 
%   [3] Starke, Gerhard. "Near-circularity for the rational Zolotarev 
%   problem in the complex plane." J. of approx. theory 70.1 (1992): 115-130.
%
%   [4] Townsend, Alex, and Heather Wilber. "On the singular values of matrices 
%   with high displacement rank." arXiv preprint arXiv:1712.05864 (2017).
%













