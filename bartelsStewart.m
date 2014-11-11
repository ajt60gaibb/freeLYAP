function X = bartelsStewart(A, B, C, D, E, xSplit, ySplit, tol)
%BARTELSSTEWART   Solution to generalized Sylvester matrix equation. 
% 
% Computes the solution to the Sylvester equation
%
%         AXB^T + CXD^T = E
%
% by using the Bartels--Stewart algorithm, see 
%
% J. D. Gardiner, A. J. Laub, J. J. Amato, & C. B. Moler, Solution of the
% Sylvester matrix equation AXB^T + CXD^T = E, ACM Transactions on Mathematical
% Software (TOMS), 18(2), 223-231.
%
% Note that if B or C are empty, they are assumed to be identity mtrices of the
% appropriate size. This is more efficient that passing identity matrices.
% 
% This Bartels--Stewart solver also takes information xsplit, ysplit so that if
% possible it decouples the even and odd modes.

% Copyright 2014 by The University of Oxford and The Chebfun2 Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Parse inputs:
if ( nargin < 6 )
    xSplit = false;
end
if ( nargin < 7 )
    xSplit = false;
end
if ( nargin < 8 )
    % Fixed tolerance:
    tol = 10*eps; 
end

% If the RHS is zero then the solution is the zero solution (assuming
% uniqueness).
if ( norm(E) < tol )
    X = zeros(size(E));
    return
end

if ( isempty(B) && isempty(C) )
%     X = sylv(A, D, E);
    X = lyap(A, D, E);
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
    [Z1, P] = schur(A, 'complex');
    Q1 = Z1';
%     S = speye(n, m);
    S = eye(n, m);
elseif ( ySplit )
    % This is equivalent to qz(A,C), but faster.
    [P, S, Q1, Z1] = qzsplit(A, C); 
else
    [P, S, Q1, Z1] = qz(A, C);  
end
% We enforce P and S as upper triangular because they should be (up to rounding
% errors) and we need to do back substitution with them.
P = triu(P); 
S = triu(S);

% If the PDE is even/odd in the y-direction then we can split (further) into
% double as many subproblems.
if ( isempty(B) )
    [Z2, T] = schur(D, 'complex');
    Q2 = Z2';
%     R = speye(m, n);
    R = eye(m, n);
elseif ( xSplit )
    % Faster QZ when even/odd modes decouple in x-direction: 
    [T, R, Q2, Z2] = qzsplit(D, B);
else
    % QZ does not take matrices in sparse format:
    [T, R, Q2, Z2] = qz(D, B);
end

% Now use the generalised Bartels--Stewart solver found in Gardiner et al.
% (1992). The Sylvester matrix equation now contains quasi upper-triangular
% matrices and we can do a backwards substitution.

% transform the righthand side.
F = Q1*E*Q2.';

% Do a backwards substitution type algorithm to construct the solution.
k = n;
PY = zeros(m);
SY = zeros(m);

% Construct columns n,n-1,...,3,2 of the transformed solution.  The first
% column is treated as special at the end.
while ( k > 1 )
    
    % There are two cases, either the subdiagonal contains a zero, i.e.,
    % T(k,k-1) = 0 and then it is simply a backwards substitution, or T(k,k-1)
    % ~= 0 and we solve a 2x2 system.
    
    if ( T(k,k-1) == 0 )
        % Simple case (Usually end up here).
        
        if ( k < n )    
            PY(:,k+1) = P*Y(:,k+1);
            SY(:,k+1) = S*Y(:,k+1);    
            jj = (k+1):n;
            rhs = F(:,k) - PY(:,jj)*R(k,jj).' - SY(:,jj)*T(k,jj).';
        else
            rhs = F(:,k);
        end 
        
        % Find the kth column of the transformed solution.
        tmp = (P + (T(k,k)/R(k,k))*S); % <- Divide both sides by R_kk for speed.
        rhs = rhs/R(k,k);
        Y(:,k) = tmp \ rhs;
        
        % Go to next column
        k = k - 1;
        
    else

        % This is a straight copy from the Gardiner et al. paper, and just
        % solves for two columns at once. (Works because of quasi-triangular
        % matrices.)
        
        % Operator reduction.
        rhs1 = F(:,k-1);
        rhs2 = F(:,k);
        
        for jj = (k+1):n
            yj = Y(:,jj);
            Pyj = P*yj;
            Syj = S*yj;
            rhs1 = rhs1 - R(k-1,jj)*Pyj - T(k-1,jj)*Syj;
            rhs2 = rhs2 - R(k,jj)*Pyj - T(k,jj)*Syj;
        end
        
        % 2 by 2 system.
        top = 1:n;
        bot = (n+1):(2*n);
%         SM = zeros(2*n);
%         SM(top,top) = R(k-1,k-1)*P + T(k-1,k-1)*S;
%         SM(top,bot) = R(k-1,k)*P + T(k-1,k)*S;
%         SM(bot,top) = R(k,k-1)*P + T(k,k-1)*S;
%         SM(bot,bot) = R(k,k)*P + T(k,k)*S;

        SM = [R(k-1,k-1)*P + T(k-1,k-1)*S , R(k-1,k)*P + T(k-1,k)*S ;
              R(k,k-1)*P + T(k,k-1)*S     , R(k,k)*P + T(k,k)*S     ];
          
%         % Permute the columns and rows: 
%         SPermuted = zeros(2*n);
%         SPermuted(1:2:2*n,1:2:2*n) = SM(top, top); 
%         SPermuted(2:2:2*n,2:2:2*n) = SM(bot, bot);

        % Solve.
        UM = SM \ [rhs1 ; rhs2];
        
        Y(:,k-1)  = UM(top); 
        Y(:,k)    = UM(bot);
        PY(:,k)   = P*Y(:,k);
        PY(:,k-1) = P*Y(:,k-1);
        SY(:,k)   = S*Y(:,k); 
        SY(:,k-1) = S*Y(:,k-1);

        % We solved for two columns so go two columns further.
        k = k - 2;
        
    end
    
end

if ( k == 1 )
    % Now we have just the first column to compute.
    PY(:,2) = P*Y(:,2);
    SY(:,2) = S*Y(:,2);
    jj = 2:n;
    rhs = F(:,1)  - PY(:,jj)*R(1,jj).' - SY(:,jj)*T(1,jj).';
    Y(:,1) = (R(1,1)*P + T(1,1)*S) \ rhs;
end

% We have now computed the transformed solution so we just transform it back.
X = Z1*Y*Z2.';

end

function [P, S, Q1, Z1] = qzsplit(A, C)
%QZSPLIT   A faster qz factorisation for problems that decouple.
%
% This is equivalent to standard QZ, except we take account of symmetry to
% reduce the computational requirements.

% Do the QZ by splitting the problem into two subproblems. 

n = size(A, 1);

odd  = 1:2:n;
even = 2:2:n;

A1 = A(odd, odd); 
C1 = C(odd, odd);
[P1, S1, Q1, Z1] = qz(A1, C1);

A2 = A(even, even); 
C2 = C(even, even);
[P2, S2, Q2, Z2] = qz(A2, C2);

[P, S, Q1, Z1] = reform(P1, P2, S1, S2, Q1, Q2, Z1, Z2);

end

function [P, S, Q, Z] = reform(P1, P2, S1, S2, Q1, Q2, Z1, Z2)
%REFORM   Recombine subproblems to form the QZ factorization. 

% Determine the indedxing:
hf1 = size(P1, 1);
n   = 2*hf1 - 1;

top = 1:hf1;
bot = (hf1+1):n;
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