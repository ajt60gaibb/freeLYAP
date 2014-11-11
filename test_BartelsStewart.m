% function pass = test_BartelsStewart(  )
% Test generalized Sylvester matrix equation solver

tol = 1e-8;

%%

n = 10; 
rng(0)
A = rand(n); B = rand(n); C = rand(n); D = rand(n); X = rand(n); 

E = A * X * B.' + C * X * D.';
Y = bartelsStewart(A, B, C, D, E, 0, 0); 

err(1) = norm( Y - X );
pass(1) = err(1) < tol; 

%%

A = rand(n) + 1i*rand(n); B = rand(n) + 1i*rand(n); 
C = rand(n) + 1i*rand(n); D = rand(n) + 1i*rand(n); 
X = rand(n) + 1i*rand(n);

E = A * X * B.' + C * X * D.';
Y = bartelsStewart(A, B, C, D, E, 0, 0); 

err(2) = norm( Y - X );
pass(2) = err(2) < tol; 

%%

A = rand(n); 
B = eye(n); 
C = eye(n);
D = rand(n); 
X = rand(n);

E = A * X * B.' + C * X * D.';
Y = bartelsStewart(A, [], [], D, E, 0, 0); 

err(3) = norm( Y - X );
pass(3) = err(3) < tol; 

%%

A = rand(n) + 1i*rand(n); 
B = eye(n); 
C = eye(n);
D = rand(n) + 1i*rand(n); 
X = rand(n) + 1i*rand(n);

E = A * X * B.' + C * X * D.';
Y = bartelsStewart(A, [], [], D, E, 0, 0); 

err(4) = norm( Y - X );
pass(4) = err(4) < tol; 

%%

% err.'

% end