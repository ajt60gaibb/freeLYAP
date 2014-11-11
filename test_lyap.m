function pass = test_lyap()

% Test LYAP against the built-in (i.e., Control Toolbox) implementation.

tol = 1e-10;
n = 9;
m = n+1;

currDir = cd;
lyapDir = fileparts(which('sylvslv'));

try

    c = [0, 1, 0, 1];
    err = zeros(2, 3);
    for j = 1:2

        rng(0)
        A = rand(n) + c(j)*1i*rand(n);
        Q = rand(n) + c(j)*1i*rand(n);
        B = rand(m) + c(j)*1i*rand(m);
        C = rand(n,m) + c(j)*1i*rand(n,m);
        E = rand(n) + c(j)*1i*rand(n);

        cd(lyapDir)
        U = lyap(A, Q);
        cd(currDir)
        V = lyap(A, Q);
        err(j,1) = norm(U - V, inf);

        cd(lyapDir)
        U = lyap(A, B, C);
        cd(currDir)
        V = lyap(A, B, C);
        err(j,2) = norm(U - V, inf);

        % Built-in LYAP is fussy in this mode:
        Q = Q + Q';
        E = real(E);
        A = real(A);
        Q = real(Q);
        cd(lyapDir)
        U = lyap(A, Q, [], E);
        cd(currDir)
        V = lyap(A, Q, [], E);
        err(j,3) = norm(U - V, inf);

    end

    pass = err < tol;
    
catch ME
    
    cd(currDir)
    rethrow(ME);
    
end

end
