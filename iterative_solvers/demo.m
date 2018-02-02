%Example 5: Poisson's equation on a square. 
%
% In this example, we use finite differences and fiadi to solve Poisson's 
% equation on the square [-1 1]^2 in low rank form. 
%
% ADI-based methods can also used to develop pectrally accurate, 
% optimal complexity Poisson solvers on rectangles, cylinders, cubes, 
% and solid spheres (See []). 
%

% written by Heather Wilber, Jan. 2018


%%
% We will solve lap(u) = f with homogenous boundary conditions. 
% Using second-order central finite differences, we define
% a differentiation operator D and discretize the above equation as Lyapunov matrix 
% equation [3]: DX + XD' = F.  Here, F is a sample of f on an equispaced
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
% The right hand side is sampled from a smooth function. We plot the
% singular values of F to show that they decay rapidly. 
tol = 1e-8; 
gam = 200; 
f = @(x,y)  1/200*exp( ((x.^(gam) + y.^(gam)).^(1/gam) )); 
 F = f(xx,yy); 
 F = F(2:end-1, 2:end-1);
 [U, S, V] = svd(F);
 s = s/s(1,1); 
 s = s(s > tol); 
 figure(1)
 semilogy(1:length(s), s, 'ko')
 title('Singular value decay of F')

%%
% To solve DX + XD' = F, i.e., AX - X(-D') = F with fiadi, we need to bound 
% the eigenvalues of D and -D'. The eigenvalues of D are real, negative and 
% bounded in I1:

I1 = [-(m+1)^2 -1]; 

%The eigenvalues of -D' are therefore bounded in I2: 

I2 = [1 (m+1)^2];

%%
% To use fiadi, we compute a set of shift parameters: 

adi_tol = 1e-10; 
[P, Q] = getshifts_fiadi([I1 I2],D, -D', U, S, V, 'tol', adi_tol); 

%%
% Now we can solve for X in low rank form: 

[Z, D, Y] = fiadi(D, -D', U, S, V, P, Q, adi_tol); 

%%
% Here's what the solution looks like (boundaries omitted): 
surf(Z*D*Y'), view(2), shading flat



%we check the accuracy of fiadi using lyap 




