function  test_adi_functions()
% test the basic performance of the adi, fadi, and fiadi functions.

tol = 1e-11;
tolt = 1e-10; 

%%
% Part 1: test with real eigenvalues: 
n = 256; 
i = 1; 
  
a = -100; b = -.01; c = .01; d = 20; 
A = rand(n,n);
[Q, ~] = qr(A); 
D = linspace(a,b,n); 
A = Q*diag(D)*Q'; % A is normal with eigenvalues in [a b]

B = rand(n,n); 
[Q, ~] = qr(B); 
D = linspace(c,d,n); 
B = Q*diag(D)*Q'; % B is normal with eigenvalues in [c d]. 

r = 1; 
U = rand(n,r); 
V = rand(n,r); 
F = U*V';

%%

I = [a b c d]; 
Xtrue = lyap(A, -B, -F);
NXtrue = norm(Xtrue); 

%test get shifts works: 
     
%(1) 
[p,~] = getshifts_adi(I, 5);
pass(i) = (length(p) ==5); 

%(2)
i = i+1;
pass(i) = 0; 
try [p,q] = getshifts_adi( [-10 12 3 25], 5); 
catch
    pass(i) =1; 
end 

%(3) 
i = i+1;
pass(i) = 0; 
try [p,q] = getshifts_adi( [-10 12 3 5], 5); 
catch
    pass(i) =1; 
end

% test various adi inputs: 
[p, q] = getshifts_adi(I, 'tol', tol); 

%(4)
i = i + 1; 
X = adi(A, B, F); 
pass(i) = ( norm(X - Xtrue)/NXtrue < tolt );

%(5)
i = i + 1;
X = adi(A, B, F, tol);
pass(i) = ( norm(X - Xtrue)/NXtrue < tolt ); 

%(6)
i = i + 1;
X = adi(A, B, F, p, q); 
pass(i) = ( norm(X - Xtrue)/NXtrue < tolt );  

% test various fadi inputs

%(7)
i = i + 1; 
[Z, D, Y] = fadi(A, B, U, V); 
pass(i) = ( norm(Xtrue - Z*D*Y')/NXtrue < tolt ); 

%(8)
i = i + 1; 
[Z, D, Y] = fadi(A, B, U, V, tol); 
pass(i) = ( norm(Xtrue - Z*D*Y')/NXtrue < tolt ); 

%(9)
i = i + 1; 
[Z, D, Y] = fadi(A, B, U, V, p, q); 
pass(i) = ( norm(Xtrue - Z*D*Y')/NXtrue < tolt ); 

%test with a higher rank rhs: 
r = 10; 
U = rand(n,r); 
V = rand(n,r); 
F = U*V';
Xtrue = lyap(A, -B, -F);
NXtrue = norm(Xtrue);

%(10)
i = i + 1; 
[Z, D, Y] = fadi(A, B, U, V);
pass(i) = ( norm(Xtrue - Z*D*Y')/NXtrue < tolt ); 

%(11)
i = i + 1; 
[Z, D, Y] = fadi(A, B, U, V, p, q); %uncompressed
pass(i) = ( norm(Xtrue - Z*D*Y')/NXtrue < tolt ); 

%(12)
i = i + 1; %with compression
[Z, D, Y] = fadi(A, B, U, V, p, q, tol);
pass(i) = ( norm(Xtrue - Z*D*Y')/NXtrue < tolt ); 

%%
% test fiadi: 
r = 100; 
U = rand(n,r); 
V = rand(n,r); 
[U, S, V] = lowrank_svd(U, V); 

rho = 2; 
s = rho.^(-(1:r)); 
S = diag(s); 

Xtrue = lyap(A, -B, -U*S*V');
NXtrue = norm(Xtrue);

%%
% test getshifts

%(13)
i = i+1;
pass(i) = 0; 
try [p,q] = getshifts_fiadi( [-10 12 3 5], A, B, U, S, V, tol); 
catch
    pass(i) =1; 
end

% test various fiadi input options: 

%(14)
i = i + 1; 
[P, Q]= getshifts_fiadi(I,A,B,U,S,V, tol);
[Z, D, Y] = fiadi(A, B, U, S, V, P, Q,tol); 
pass(i) =  ( norm(Xtrue - Z*D*Y')/NXtrue < tolt ); 

%(15)
i = i+1;
[Z, D, Y] = fiadi(A, B, U, S, V, P, Q, 'nocompress');
pass(i) = ( norm(Xtrue - Z*D*Y')/NXtrue < tolt ); 

%(16)
i = i + 1; 
[Z, D, Y] = fiadi(A, B, U, S, V, tol); %compute shift parameters auto.
pass(i) = ( norm(Xtrue - Z*D*Y')/NXtrue < 5*tolt );

%(17)
i = i + 1; 
[Z, D, Y] = fiadi(A, B, U, S, V); %compute shift parameters auto.
pass(i) = ( norm(Xtrue - Z*D*Y')/NXtrue < 5*tolt );


% a few more tests on adi and fadi: 
%
% mix up order of eigenvalue intervals; make sure shift parameters 
% still compute correctly: 

temp = A;  
A = B; 
B = temp; 
 
U = rand(n,1); 
V = rand(n,1); 
F = U*V';

Xtrue = lyap(A, -B, -F);
NXtrue = norm(Xtrue);
I = [ c d a b]; 
[p, q] = getshifts_adi(I,'tol', tol); 

%(18)
i = i + 1;
X = adi(A, B, F, p, q); 
pass(i) = ( norm(X - Xtrue)/NXtrue < tolt );  

%(19)
i = i + 1; 
[Z, D, Y] = fadi(A, B, U, V, p, q, tol); 
pass(i) = ( norm(Xtrue - Z*D*Y')/NXtrue < tolt ); 

%test close intervals: 
a = -100; b = -.00001; c = .00001; d = 20; 
A = rand(n,n);
[Q, ~] = qr(A); 
D = linspace(a,b,n); 
A = Q*diag(D)*Q'; 

B = rand(n,n); 
[Q, ~] = qr(B); 
D = linspace(c,d,n); 
B = Q*diag(D)*Q';

Xtrue = lyap(A, -B, -F);
NXtrue = norm(Xtrue);
I = [ a b c d]; 
[p, q] = getshifts_adi(I,'tol', tol); 

%(20)
i = i + 1;
X = adi(A, B, F, p, q); 
pass(i) = ( norm(X - Xtrue)/NXtrue < 15*tolt );  

%(21)
i = i + 1; 
[Z, D, Y] = fadi(A, B, U, V, p, q, tol); 
pass(i) = ( norm(Xtrue - Z*D*Y')/NXtrue < 10*tolt ); 

% test fiadi with close intervals: 
r = 100; 
U = rand(n,r); 
V = rand(n,r); 
[U, S, V] = lowrank_svd(U, V); 

rho = 2; 
s = rho.^(-(1:r)); 
S = diag(s); 

Xtrue = lyap(A, -B, -U*S*V');
NXtrue = norm(Xtrue);

%(22)
i = i + 1; 
[P, Q]= getshifts_fiadi(I,A,B,U,S,V, tol);
[Z, D, Y] = fiadi(A, B, U, S, V, P, Q,tol); 
pass(i) = ( norm(Xtrue - Z*D*Y')/NXtrue < 10*tolt ); 

%%
% test with complex eigenvalues: 

% set up: 
o1 = -10+0i; 
r1 = 5; 
o2 = 2 + 5i; 
r2 = 3; 

A = rand(n, n); 
[Q, ~] = qr(A,0); 
t = -pi + rand(n,1)*2*pi; 
rad = rand(n,1);
eigsA = rad*(r1).*exp(1i*t)+o1;
A = Q*diag(eigsA)*Q'; 

B = rand(n, n);
[Q, ~] = qr(B, 0); 
t = -pi + rand(n,1)*2*pi; 
rad = rand(n,1);
eigsB = rad*(r2).*exp(1i*t)+o2;
B = Q*diag(eigsB)*Q'; 

r = 1; 
U = rand(n,r); 
V = rand(n,r); 
F = U*V';

I = [o1 r1 o2 r2]; 
Xtrue = lyap(A, -B, -F);
NXtrue = norm(Xtrue); 

% test getshifts_smith: 

%(23)
i = i+1;
pass(i) = 0;
try [p,q] = getshifts_smith([1 3 2 2], 'tol', tol);
catch
    pass(i) = 1;
end

[p, q] = getshifts_smith(I, 'tol', tol); 

%(24)
i = i + 1; 
pass(i)=0; 
try X = adi(A, B, F);
catch
    pass(i) = 1;
end

%(25)
i = i + 1; 
X = adi(A, B, F, p, q); 
pass(i) = ( norm(Xtrue - X)/NXtrue < tolt ); 

%(26)
i = i + 1; 
pass(i)=0; 
try [Z, D, Y] = fadi(A, B, U,V);
catch
    pass(i) = 1;
end

%(27)
i = i + 1; 
[Z, D, Y] = fadi(A, B, U,V, p, q); 
pass(i) = ( norm(Xtrue - Z*D*Y')/NXtrue < tolt ); 

% test fiadi: 
r = 100; 
U = rand(n,r); 
V = rand(n,r); 
[U, S, V] = lowrank_svd(U, V); 

rho = 2.4; 
s = rho.^(-(1:r)); 
S = diag(s); 

Xtrue = lyap(A, -B, -U*S*V');
NXtrue = norm(Xtrue);

%(28)
i = i+1;
pass(i) = 0;
try [P,Q] = getshifts_fismith([1 3 2 2],A, B, U, S, V, tol);
catch
    pass(i) = 1;
end

%(29)
i = i+1;
pass(i) = 0;
try [Z, D, Y] = fiadi(A, B, U,S, V,tol); 
catch
    pass(i) = 1;
end

%(30)

[P,Q] = getshifts_fismith(I,A, B, U, S, V, tol);

i = i + 1; 
[Z, D, Y] = fiadi(A, B, U,S, V, P, Q, tol); 
pass(i) = ( norm(Xtrue - Z*D*Y')/NXtrue < 5*tolt ); 

%%
% summarize

if any(pass==0)
    idx = find(pass-1);
    fprintf('The following tests failed:') 
    disp(idx)
else
    fprintf('All tests passed.')

end
end






