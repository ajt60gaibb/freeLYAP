function X = lyap(A, B, C, E)
%LYAP  Solve continuous-time Lyapunov equations.
%   LYAP() aims to mirror the functionality and syntax of the LYAP() function in
%   the MATLAB Control Toolbox. It is written entirely in MATLAB and so is a
%   little slower than the implementation in the Control Toobox (which is
%   essentially a wrapper for LAPACK/SLICOT routines). However it is a little
%   more flexible.
%
%   X = LYAP(A,Q) solves the Lyapunov matrix equation:
%
%       A*X + X*A' + Q = 0
%
%   X = LYAP(A,B,C) solves the Sylvester equation:
%
%       A*X + X*B + C = 0
%
%   X = LYAP(A,Q,[],E) solves the generalized Lyapunov equation:
%
%       A*X*E' + E*X*A' + Q = 0.
%
%   Note that unlike the built-in MATLAB LYAP() routine there is no need for Q
%   to be symmetric in the final case. Also note that no balancing is performed.
%
%   For solving generalised Sylvestter equations fo the form 
%       A*X*B^T + C*X*D^T = E
%   see bartelsStewart.m.
%
% See also BARTELSSTEWART.

% Nick Hale, Nov 2014. (nick.p.hale@gmail.com)

% TODO: Balancing?

if ( nargin == 2 )
    % A*X + X*A' + B = 0
    X = sylv(A, [], B);
    if ( isreal(A) && isreal(B) )
        X = real(X);
    end
    
elseif ( nargin == 3 )
    % A*X + X*B + C = 0
    X = sylv(A, B, C);
    if ( isreal(A) && isreal(B) && isreal(C) )
        X = real(X);
    end
    
else
    % A*X*E' + E*X*A' + B = 0
    X = bartelsStewart(A, E, E, A, -B, false, false);
    if ( isreal(A) && isreal(B) && isreal(E) )
        X = real(X);
    end
    
end

end

function X = sylv(A, B, C)
%SYLV  Solve Sylvester matrix equation.
%   SYLV(A, B, C) solves the Sylvester matrix equation
%       A*X + X*B + C = 0.
%   SYLV(A, [], C) solves the Lyapunov matrix equation
%       A*X + X*A' + C = 0.

% Nick Hale, Nov 2014. (nick.p.hale@gmail.com)

% Get sizes:
m = size(A, 1); 
n = size(B, 1); 

% Compute Schur factorizations. (P and T will be upper triangular.)
[Z1, P] = schur(A, 'complex');
if ( ~isempty(B) )
    [Z2, T] = schur(B.', 'complex');
else
    % B = A';
    % [Z2, T] = schur(B.', 'complex');
    Z2 = conj(Z1);
    T = conj(P);
    n = size(A, 2);
end

% Transform the righthand side.
F = Z1'*C*conj(Z2); % <-- Not sure why conj() is needed here.

% Initialise identity matrix and storage for transformed solution.
I = speye(size(P));
Y = zeros(m, n);

% Do a backwards substitution to construct the transformed solution.
for k = n:-1:1
    jj = (k+1):n;
    rhs = F(:,k) + Y(:,jj)*T(k,jj).';
    % Find the kth column of the transformed solution.
    Y(:,k) = (P + T(k,k)*I) \ (-rhs);
end

% Transform solution back:
X = Z1*Y*Z2.';

end