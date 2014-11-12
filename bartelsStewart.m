function X = bartelsStewart(A, B, C, D, E, split)
%BARTELSSTEWART   Solution to generalized Sylvester matrix equation. 
% 
% BARTELSSTEWART(A, B, C, D, E) computes the solution to the Sylvester equation
%
%         AXB^T + CXD^T = E
%
% by using the Bartels--Stewart algorithm [1].
%
% If B or C are empty, they are assumed to be identity matrices of the
% appropriate size. This is more efficient that passing identity matrices.
% 
% BARTELSSTEWART(A, B, C, D, E, SPLIT, YSPLIT) also takes information XSPLIT,
% YSPLIT which state whether or not the problem may be decoupled into even and
% odd modes for efficiency.
%
% References:
%  [1] J. D. Gardiner, A. J. Laub, J. J. Amato, & C. B. Moler, Solution of the
%  Sylvester matrix equation AXB^T + CXD^T = E, ACM Transactions on Mathematical
%  Software (TOMS), 18(2), 223-231, 1992.

% Copyright 2014 by The University of Oxford and The Chebfun2 Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Written by Alex Townsend, ??? 2014
% Modified by Nick Hale, Nov 2014. (nick.p.hale@gmail.com)

% Parse inputs:
if ( nargin < 6 )
    split = [false ; false];
end

% Fixed tolerance:
tol = 10*eps; 

% If the RHS is zero then the solution is the zero solution (assuming
% uniqueness).
if ( norm(E) < tol )
    X = zeros(size(E));
    return
end

% Treat the standard Lyapunov and Sylvester equations as special cases:
if ( isempty(B) && isempty(C) )
    X = lyap(A, D.', -E);
    return
end

% Matrices must be sparse for QZ():
A = full(A); B = full(B); C = full(C); D = full(D);

% Solution will be a m by n matrix.
m = size(A, 1); 
n = size(D, 2); 
Y = zeros(m, n);

% If the equation is even/odd in the x-direction then we can split the problem
% into two subproblems.
if ( isempty(C) )
%     [Z1, P] = schur(A, 'complex');
    [Z1, P] = schur(A);
    Q1 = Z1';
    S = eye(n, m);
elseif ( split(2) )
    % This is equivalent to qz(A,C), but faster.
    [P, S, Q1, Z1] = qzSplit(A, C); 
else
    [P, S, Q1, Z1] = qz(A, C);  
end

% We enforce P and S as upper triangular because they should be (up to rounding
% errors) and we need to do back substitution with them.
% P = triu(P); 
% S = triu(S);

% If the PDE is even/odd in the y-direction then we can split (further) into
% double as many subproblems.
if ( isempty(B) )
%     [Z2, T] = schur(D, 'complex');
    [Z2, T] = schur(D);
    Q2 = Z2';
    R = eye(m, n);
elseif ( split(1) )
    % Faster QZ when even/odd modes decouple in x-direction: 
    [T, R, Q2, Z2] = qzSplit(D, B);
else
    % QZ does not take matrices in sparse format:
    [T, R, Q2, Z2] = qz(D, B);
end

% Now use the generalised Bartels--Stewart solver found in Gardiner et al.
% (1992). The Sylvester matrix equation now contains quasi upper-triangular
% matrices and we can do a backwards substitution.

% transform the righthand side.
F = Q1*E*Q2.';

% Initialise S*Y and P*Y factors:
PY = zeros(m, n);
SY = zeros(m, n);

% Do a backwards substitution type algorithm to construct the solution.
k = n;

% Construct columns n,n-1,...,3,2,1 of the transformed solution.
while ( k > 0 )
    
    % There are two cases, either the subdiagonal contains a zero, i.e.,
    % T(k,k-1) = 0 and then it is simply a backwards substitution, or T(k,k-1)
    % ~= 0 and we solve a 2mx2m system.
    
    if ( k == 1 || T(k,k-1) == 0 )
        % Simple case (Usually end up here).
        
        jj = (k+1):n;
        rhs = F(:,k) - PY(:,jj)*R(k,jj).' - SY(:,jj)*T(k,jj).';
        
        % Find the kth column of the transformed solution.
        tmp = (P + (T(k,k)/R(k,k))*S); % <- Divide both sides by R_kk for speed.
        rhs = rhs/R(k,k);
        Y(:,k) = tmp \ rhs;
        
        % Store S*Y and P*Y factors:
        PY(:,k) = P*Y(:,k);
        SY(:,k) = S*Y(:,k);
        
        % Go to next column:
        k = k - 1;
        
    else

        % This is a straight copy from the Gardiner et al. paper, and just
        % solves for two columns at once. (Works because of quasi-triangular
        % matrices.)
        
        % Operator reduction.
        jj = (k+1):n;
        rhs1 = F(:,k-1) - PY(:,jj)*R(k-1,jj).' - SY(:,jj)*T(k-1,jj).';
        rhs2 = F(:,k)   - PY(:,jj)*R(k,jj).'   - SY(:,jj)*T(k,jj).';

        % 2 by 2 system.
        SM = [R(k-1,k-1)*P + T(k-1,k-1)*S , R(k-1,k)*P + T(k-1,k)*S ;
              T(k,k-1)*S                  , R(k,k)*P + T(k,k)*S     ];
        
        % Solve.
%         UM = SM \ [rhs1 ; rhs2];

        % Solve (permute the columns and rows):
        idx = reshape([(1:n) ; (n+1:2*n)], 2*n, 1);
        rhs = [rhs1 ; rhs2];
        UM = SM(idx,idx) \ rhs(idx);
        UM(idx) = UM;
        
        % Store S*Y and P*Y factors:
        Y(:,k-1:k) = reshape(UM, n, 2);
        PY(:,k-1:k) = P*Y(:,k-1:k);
        SY(:,k-1:k) = S*Y(:,k-1:k);

        % We solved for two columns so go two columns farther.
        k = k - 2;
        
    end
    
end

% We have now computed the transformed solution so we just transform it back.
X = Z1*Y*Z2.';

end

function [P, S, Q1, Z1] = qzSplit(A, C)
%QZSPLIT   A faster QZ factorisation for problems that decouple.
%   QZSPLIT() is equivalent to standard QZ, except we take account of symmetry
%   to reduce the computational requirements.

% Matrix size (square):
n = size(A, 1);

% Do the QZ by splitting the problem into two subproblems. 

% Odd part:
odd  = 1:2:n;
A1 = A(odd, odd); 
C1 = C(odd, odd);
[P1, S1, Q1, Z1] = qz(A1, C1);

% Even part:
even = 2:2:n;
A2 = A(even, even); 
C2 = C(even, even);
[P2, S2, Q2, Z2] = qz(A2, C2);

% Recombine:
[P, S, Q1, Z1] = qzRecombine(P1, P2, S1, S2, Q1, Q2, Z1, Z2);

end

function [P, S, Q, Z] = qzRecombine(P1, P2, S1, S2, Q1, Q2, Z1, Z2)
%QZRECOMBINE  Recombine subproblems to form the QZ factorization. 

% Determine the indexing:
hf1  = size(P1, 1);
n    = 2*hf1 - 1;
top  = 1:hf1;
bot  = (hf1+1):n;
odd  = 1:2:n;
even = 2:2:n;

% Push the subproblem back together:
P = blkdiag(P1, P2);
% P = zeros(n);
% P(top,top) = P1; 
% P(bot,bot) = P2;

S = blkdiag(S1, S2);
% S = zeros(n);
% S(top,top) = S1; 
% S(bot,bot) = S2;

Q = zeros(n);
Q(top,odd) = Q1; 
Q(bot,even) = Q2;

Z = zeros(n);
Z(odd, top) = Z1; 
Z(even,bot) = Z2;

end