clc
% Test LYAP against the built-in (i.e., Control Toolbox) implementation.

n = 100;
m = n+1;
loop = 10;

currDir = cd;
lyapDir = fileparts(which('sylvslv'));

rng(0)
A = rand(n);
Q = rand(n);
B = rand(m);
C = rand(n,m);
E = rand(n);

cd(lyapDir)
fprintf('built-in real AQ:\n')
tic
for k = 1:loop
    U = lyap(A, Q);
end
toc
%%
cd(currDir)
fprintf('new real AQ:\n')
tic
for k = 1:loop
    U = lyap(A, Q);
end
toc

disp(' ')

%%

cd(lyapDir)
fprintf('built-in real ABC:\n')
tic
for k = 1:loop
    U = lyap(A, B, C);
end
toc
cd(currDir)
fprintf('new real ABC:\n')
tic
for k = 1:loop
    U = lyap(A, B, C);
end
toc

disp(' ')

%%
% Built-in LYAP is fussy in this mode:
Q = Q + Q';
cd(lyapDir)
fprintf('built-in real AQ[]E:\n')
tic
for k = 1:loop
    U = lyap(A, Q, [], E);
end
toc
cd(currDir)
fprintf('new real AQ[]E:\n')
tic
for k = 1:loop
    U = lyap(A, Q, [], E);
end
toc

disp(' ')

%%

rng(0)
A = rand(n) + 1i*rand(n); Q = rand(n) + 1i*rand(n);
B = rand(m) + 1i*rand(m); C = rand(n,m) + 1i*rand(n,m);
E = rand(n) + 1i*rand(n);

cd(lyapDir)
fprintf('built-in complex AQ:\n')
tic
for k = 1:loop
    U = lyap(A, Q);
end
toc
cd(currDir)
fprintf('new complex AQ:\n')
tic
for k = 1:loop
    U = lyap(A, Q);
end
toc

disp(' ')

%%

cd(lyapDir)
fprintf('built-in complex ABC:\n')
tic
for k = 1:loop
    U = lyap(A, B, C);
end
toc
cd(currDir)
fprintf('new complex ABC:\n')
tic
for k = 1:loop
    U = lyap(A, B, C);
end
toc
