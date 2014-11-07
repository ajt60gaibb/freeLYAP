function pass = test_lyap()

% Test LYAP against the built-in (i.e., Control Toolbox) implementation.

tol = 1e-10;
n = 10;

currDir = cd;
lyapDir = fileparts(which('sylvslv'));

try

    c = [0, 1];
    err = zeros(2, 3);
    for j = 1:2

        seedRNG(0)
        A = rand(n) + c(j)*1i*rand(n);
        B = rand(n) + c(j)*1i*rand(n);
        C = rand(n) + c(j)*1i*rand(n);
        E = rand(n) + c(j)*1i*rand(n);

        cd(lyapDir)
        U = lyap(A, B);
        cd(currDir)
        V = lyap(A, B);
        err(j,1) = norm(U - V, inf);

        cd(lyapDir)
        U = lyap(A, B, C);
        cd(currDir)
        V = lyap(A, B, C);
        err(j,2) = norm(U - V, inf);

        % Built-in LYAP is fussy in this mode:
        B = B + B';
        E = real(E);
        A = real(A);
        B = real(B);
        cd(lyapDir)
        U = lyap(A, B, [], E);
        cd(currDir)
        V = lyap(A, B, [], E);
        err(j,3) = norm(U - V, inf);

    end

    pass = err < tol;
    
catch ME
    
    cd(currDir)
    rethrow(ME);
    
end

end
