%Example 5: Poisson's equation on a square. 
%
% In this example, we use finite differences and fiadi to solve Poisson's 
% equation on the square [-1 1]^2. We show how to find the solution in low 
% rank form (with fiadi), and also explicitly (with adi). 
%
% ADI-based methods can also used to develop spectrally accurate, 
% optimal complexity Poisson solvers on rectangles, cylinders, cubes, 
% and solid spheres (See [1,2]). 
%

% written by Heather Wilber, Jan. 2018


%%
% We will solve lap(u) = f with homogenous boundary conditions. 
% Using second-order central finite differences, we define
% a differentiation operator D and discretize the above equation as Lyapunov matrix 
% equation [1]: DX + XD' = F.  Here, F is a sample of f on an equispaced
% grid.
% 

%%
% Set up the grid: 
clear all
m = 512; 
a = -1; 
b = 1;

h = (b-a)/(m+1); %grid spacing parameter. 
[xx,yy] = meshgrid(linspace(a,b,m+2),linspace(a,b,m+2));

%%
% set up the differentiation operator:
D = 1/h^2*spdiags([1*ones(m,1) -2*ones(m,1) 1*ones(m,1)],[-1 0 1], m,m);

%%
% The right hand side is sampled from a smooth function and
% truncated. We plot the singular values of F to show that they decay rapidly. 
tol = 1e-8; 
gam = 200; 
f = @(x,y)  1/200*exp( ((x.^(gam) + y.^(gam)).^(1/gam) )); 
 F = f(xx,yy); 
 F = F(2:end-1, 2:end-1);
 [U, S, V] = svd(F);
 s = diag(S); s = s/s(1,1); s = s(s > tol); 
 S = diag(s); U = U(:,1:length(s)); V = V(:,1:length(s)); 
 F = U*S*V';
 semilogy(1:length(s), s, 'ko'), title('Singular value decay of F')
 set(gca,'fontsize',18)

%%
% To solve DX + XD' = USV', i.e., AX - X(-D') = USV' with fiadi, we need to bound 
% the eigenvalues of D and -D'. The eigenvalues of D are real, negative and 
% bounded in I1 (see [1]):

I1 = [-(m+1)^2 -1]; 

%The eigenvalues of -D' are therefore bounded in I2: 

I2 = [1 (m+1)^2];

%%
% To use fiadi, we compute a set of shift parameters: 

adi_tol = 1e-8; 
[P, Q] = getshifts_fiadi([I1 I2],D, -D', U, S, V, adi_tol); 

%%
% Now we can solve for X in low rank form: 

[Z, DD, Y] = fiadi(D, -D', U, S, V, P, Q, adi_tol); 

%%
% Here's what the solution looks like (boundaries omitted): 
surf(Z*DD*Y'), view(2), shading flat

% The rank of the solution = # columns of Z: 

length(Z(1,:))

%%
% Here, we don't take advantage of the fact
% that D is diagonalized by the DFT of type 1. A fast, ADI-based low rank 
% Poisson solver using the DFT is discussed in [2].

% If we aren't interested in a low rank solution, we can find X explicitly
% using adi. First, we get adi shift parameters: 

[p, q] = getshifts_adi([I1 I2],'tol', adi_tol); 

%%
% now we simply call adi: 

X = adi(D, -D', F, p, q); 

%%
% The low rank and explicit solutions are approximately the same: 

norm(X - Z*DD*Y')

%% References
%
%   [1] Fortunato, Daniel, and Alex Townsend. 
%   "Fast Poisson solvers for spectral methods." arXiv preprint arXiv:1710.11259 (2017).
%
%   [2] Townsend, Alex, and Heather Wilber. 
%   "On the singular values of matrices with high displacement rank." 
%   arXiv preprint arXiv:1712.05864 (2017).
%





