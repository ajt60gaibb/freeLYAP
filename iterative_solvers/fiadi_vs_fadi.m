% a test to see when fiadi is beneficial: 

%first test: interval eigs: 
%%
% set up:
clear all
nn = [100,512]; 
a = 0.0001; 
b = 5; 
c = -100; 
d = -.000001; 
I = [a b c d];
tol = 1e-15; 
%%
% test

for j = 1:length(nn)
    n = nn(j); 
    %A = rand(n,n); 
    %[Q, ~] = qr(A); 
    dg  = linspace(a,b, n); 
    A = diag(dg); 
   

    %B = rand(n,n); 
    %[Q, ~] = qr(B); 
    dg = linspace(c,d, n); 
    B = diag(dg); % B is normal, has real eigenvalues in the interval [c d]
  
    % set up RHS
    r = n; 
    U = rand(n,r); 
    V = rand(n,r); 
    [U,~,V] = lowrank_svd(U, V); 

    %replace S with rapidly decaying sequence of singular values:
    rho = 4;  
    s = rho.^(-sqrt(1:r)); 
    S = diag(s); 
    % Now we have a RHS F = U*S*V'. 
    
%     % keep only sv > eps
idx = find(diag(S)>tol*S(1,1), 1, 'last');
U = U(:,1:idx); 
V = V(:,1:idx); 
S = S(1:idx, 1:idx); 

    % do timings: 
    
    %fadi w/compression: 
    s = tic;
    UU = U*sqrt(S); 
    VV = V*sqrt(S);  
    [p, q] = getshifts_adi(I, 'tol', tol); 
    [Z, D, Y] = fadi(A, B, UU, VV, p,q, tol); 
    time_fadi1(j) = toc(s); 
    
    %fiadi: 
    s = tic;
    [P, Q] = getshifts_fiadi(I, A, B, U, S, V, tol); 
    [Z, D, Y] = fiadi(A, B, U, S, V, P, Q,tol);   
    time_fiadi(j) = toc(s); 
    j
end
    
%%
% plots


    
    


