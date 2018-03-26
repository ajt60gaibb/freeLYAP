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
%   For solving generalised Sylvestter equations of the form 
%       A*X*B^T + C*X*D^T = E
%   see bartelsStewart.m.
%
%   When solving for multiple righthand sides or with various constant diagonal
%   shifts in A and/or B, it can be beneficial to precompute the Schur
%   factorizations. In particular, V = lyap(A, C) can be computed as
%
%       [U, T] = schur(A, 'complex');
%       V = U * lyap(T, U'*C*U) * U';
%
%   and V = lyap(A, B, C) as
%
%       [UA, TA] = schur(A, 'complex');
%       [UB, TB] = schur(B, 'complex');
%       V = UA * lyap(TA, TB, UA'*C*UB) * UB';
%
% See also SCHUR, BARTELSSTEWART.

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
%     X = bartelsStewart(A, E, E, A, -B);
    X = bartelsStewart(A, E, [], [], -B);
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
[m, n] = size(C); 

% Compute Schur factorizations. TA will be upper triangular. TB will be upper or
% lower. If TB is upper triangular then we backward solve; if it's lower
% triangular then do forward solve.
[ZA, TA] = schur(A, 'complex');
if ( isempty(B) || isequal(A, B') ) % <-- Should we check for isequal here?
    ZB = ZA;                        % (Since user _should_ have passed B = [].)
    TB = TA';
    n = size(A, 2);
    solve_direction = 'backward';
elseif ( isequal(A, B.') )          % <-- We _should_ here as no way to specify.
    ZB = conj(ZA);
    TB = TA.';
    solve_direction = 'backward';
else
    % We must also compute the schur factorization of B and forward solve;
    [ZB, TB] = schur(B, 'complex');
    solve_direction = 'forward';
end

% Transform the righthand side.
F = ZA' * C * ZB;

% Symmetric case is trivial!
if ( isdiag(TA) && isdiag(TB) )
    L = -1./(diag(TA) + diag(TB).');
    X = ZA*(L.*F)*ZB';
    return 
end

% Initialise storage for transformed solution.
Y = zeros(m, n);

% Diagonal mask (for speed in shifted solves):
idx = diag(true(m, 1));
p = diag(TA);

% Forwards or backwards?
if ( strcmp(solve_direction, 'backward') )
    kk = n:-1:1;
else
    kk = 1:n;
end

% Do a backwards/forwards substitution to construct the transformed solution.
for k = kk
    rhs = F(:,k) + Y*TB(:,k); % <-- More efficient multiplication.
    % Find the kth column of the transformed solution.
    TA(idx) = p + TB(k,k);    % <-- Diagonal shift. More efficient this way.
    Y(:,k) = TA \ (-rhs);
end

% Transform solution back:
X = ZA * Y * ZB';

end