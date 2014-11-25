function pass = test_matrixEquationSolver() 
% Test matrixEquationSolver.m 

tol = 1e-8; 

n = 10; 
m = 12; 

%% Test 1

A = randn( n ) + 1i*randn(n); 
B = randn( m ) + 1i*randn(m); 
EXACT = randn( n, m); 
F = A * EXACT * B.'; 

X = matrixEquationSolver({A}, {B}, F);
pass(1) = ( norm( X - EXACT ) < tol ); 

%% Test 2

C = randn( n ) + 1i*randn(n); 
D = randn( m ) + 1i*randn(m); 
EXACT = randn( n, m); 
F = A * EXACT * B.' + C * EXACT * D.'; 
X = matrixEquationSolver({A,C}, {B,D}, F);
pass(2) = ( norm( X - EXACT ) < tol ); 

%% Test 3

E = randn( n ) + 1i*randn(n) ; 
G = randn( m ) + 1i*randn(m); 
EXACT = randn( n, m); 
F = A * EXACT * B.' + C * EXACT * D.' + E * EXACT * G.'; 
X = matrixEquationSolver({A,C,E}, {B,D,G}, F);
pass(3) = ( norm( X - EXACT ) < tol ); 

%% Test 4

K = 10;  % number of terms: 
EXACT = randn( n, m ); 
F = zeros( n, m ); 
A = cell(K,1); B = cell(K,1); 
for j = 1:K 
    A{j} = randn( n ) + 1i*randn( n ); 
    B{j} = randn( m ) + 1i*randn( m ); 
    F = F + A{j} * EXACT * B{j}.'; 
end
X = matrixEquationSolver(A, B, F);
pass(4) = ( norm( X - EXACT ) < tol ); 

end