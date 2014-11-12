function pass = test_BartelsStewart(  )
% Test generalized Sylvester matrix equation solver

tol = 1e-8;

%%

n = 10; m = n-1;
rng(0)
A = rand(m); B = rand(n); C = rand(m); D = rand(n); X = rand(m,n); 

E = A * X * B.' + C * X * D.';
Y = bartelsStewart(A, B, C, D, E); 

err(1) = norm( Y - X );
pass(1) = err(1) < tol; 

%%

A = rand(m) + 1i*rand(m); B = rand(n) + 1i*rand(n); 
C = rand(m) + 1i*rand(m); D = rand(n) + 1i*rand(n); 
X = rand(m,n) + 1i*rand(m,n);

E = A * X * B.' + C * X * D.';
Y = bartelsStewart(A, B, C, D, E); 

err(2) = norm( Y - X );
pass(2) = err(2) < tol; 

%%

A = rand(m); 
B = eye(n); 
C = eye(m);
D = rand(n); 
X = rand(m,n);

E = A * X * B.' + C * X * D.';
Y = bartelsStewart(A, [], [], D, E); 

err(3) = norm( Y - X );
pass(3) = err(3) < tol; 

%%

A = rand(m) + 1i*rand(m); 
B = eye(n); 
C = eye(m);
D = rand(n) + 1i*rand(n); 
X = rand(m,n) + 1i*rand(m,n);

E = A * X * B.' + C * X * D.';
Y = bartelsStewart(A, [], [], D, E); 

err(4) = norm( Y - X );
pass(4) = err(4) < tol; 

%%

% err.'

end