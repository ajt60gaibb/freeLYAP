% These basic examples show how fadi (Example 3) and 
% fiadi (Example 4) can be used to find a low rank approximation to X in
% AX - XB = U*V' when the eigenvalues of A and B are complex. 
% This requires the use of getshifts_smith and getshifts_fismith.
% The comments discuss in more detail when it is appropriate to use these methods. 
% For a more detailed introduction to fadi and fiadi, see 
% Examples_lowrankadi12.m. 
%

% Heather Wilber, Jan. 2018

%%
% Example 3: Using fadi with complex eigenvalues
%
% We'll start by creating A, B, and a right-hand side F = U*V'.
% fadi is highly effective when rank(F) is small, 
% A and B are normal matrices, and the eigenvalues of A and B are are contained in 
% disjoint, real intervals, or disjoint disks in the complex plane (see [2]). 
% In this example, we set up A and B with eigenvalues contained in disjoint 
% disks. See Examples_lowrankadi12.m for an example with intervals. 


%%
% we choose A to have eigenvalues in a disk with center o1  and
% radius r1, and B to have eigenvalues in a disk with center o2, radius r2.

clear all
n = 100; 

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

%a plot of the eigenvalues
plot(r1*exp(1i*linspace(-pi, pi, n)) + o1, '-b', 'Linewidth', 2.5), hold on
plot(r2*exp(1i*linspace(-pi, pi, n)) + o2, '-r', 'Linewidth', 2.5)
plot(eigsA, '.b', 'Markersize', 20)
plot(eigsB, '.r', 'Markersize', 20)
legend('A', 'B', 'Location', 'Southeast')
axis equal, set(gca,'fontsize',18), xlabel('Real'),ylabel('Imaginary')
hold off

%%
% fadi works best when the right-hand side F = U*V' has a very low rank: 

rnk = 1; 
U = rand(n,rnk)*rand(rnk, n); 
V = rand(n,rnk)*rand(rnk,n); 

%%
% Using fadi to approximate X in AX-XB = UV' requires selecting a set of 
% ADI shift parameters. When the eigenvalues of A and B are 
% contained in disjoint disks, it is optimal to use a single ADI shift
% parameter repeatedly [1]. In this setting, ADI is equivalent to Smith's
% method. We can use the following command to get the
% optimal Smith's (aka ADI) shift parameter (see [1,2]): 

L = [o1 r1 o2 r2];
[p,q] = getshifts_smith(L);

%%
% To apply N fadi iterations, we use the following: 

N = 5; 
[P,Q] = getshifts_smith(L,N); 

% Here, P = p*ones(1,N), Q = q*ones(1,N). This formats the Smith's shift 
% parameter for use in fadi. The solution to AX - XB = UV' is returned 
% as X \approx Z*D*Y'. 

[Z, D, Y] = fadi(A, B, U, V, P, Q); 

% We'll check the accuracy by comparing Z*D*Y' to the solution computed via
% lyap: 

Xt = lyap(A, -B,-U*V'); 
norm(Xt - Z*D*Y')/norm(Xt)

%%
% We can also specify a tolerance parameter when computing our shifts
% that guarantees the relative accuracy of the solution: 

tol = 1e-10; 
[P,Q] = getshifts_smith(L, 'tol', tol); 
[Z, D, Y] = fadi(A, B, U, V, P, Q); 
 
norm(Xt - Z*D*Y')/norm(Xt)

%%
% Example 4: Using getshifts_fismith: 
%
% fiadi is a generalization of fadi that can be used when F in AX-XB =F
% is of moderate to high rank, but has singular values that decay rapidly. 
% We now set up a right hand side F where fiadi is useful. See
% Example_lowrankadi12.m for a more detailed discussion of fiadi. 

r = 50; 
U = rand(n,r); 
V = rand(n,r); 
[U,~,V] = lowrank_svd(U, V); 

%replace S with rapidly decaying sequence of singular values:
rho = 1.5;  
s = rho.^(-(1:r)); 
S = diag(s); 
% Now we have a RHS F = U*S*V'. 

%%
% When A and B are as in the previous example, we find the shift parameters
% to use in fiadi with the getshifts_fismith command. This command always
% requires a tolerance parameter as input: 

[P,Q] = getshifts_fismith(L, A, B, U, S, V, tol); 


% Now P and Q are each matrices, with the same value (the Smith shift
% parameter) given for each nonzero entry. This is the correct format for
% shift parameter sets passed into fiadi (see Example 2 in Examples_real_eigs.m, 
% as well as [2]).

%%
% we can now use fiadi with P and Q as our collection of shift parameters: 

[ZZ, DD, YY] = fiadi(A, B, U, S, V, P, Q, tol);


Xt = lyap(A, -B,-U*S*V'); 
norm(Xt - ZZ*DD*YY')/norm(Xt)


%%
% References: 
%
% [1] Starke, Gerhard. "Near-circularity for the rational Zolotarev 
% problem in the complex plane." J. of approx. theory 70.1 (1992): 115-130.
%
% [2] Townsend, Alex, and Heather Wilber. "On the singular values of matrices 
% with high displacement rank." arXiv preprint arXiv:1712.05864 (2017).
%








