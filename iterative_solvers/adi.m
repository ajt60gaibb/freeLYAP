function X = adi(A,B,F,varargin)
%
% adi(A, B, F, p, q)
% approximately solves AX - XB = F using ADI with ADI shift parameters 
% provided by vectors p, q.
%
% adi(A, B, F) If A and B have real eigenvalues that are contained in 
% disjoint intervals, the optimal shift parameters are computed automatically
% and the problem is solved to a relative accuracy of approximately machine 
% precision. 
%
% adi(A, B, F, tol) is as above, except the relative accuracy of the 
% of the solution is specified by tol. 
%
% See getshifts_adi and getshifts_smith for help computing shift parameters. 
%
% References: 
%
% [1] Lu, An, and Eugene L. Wachspress. 
%  "Solution of Lyapunov equations by alternating direction implicit iteration." 
%  Comp. & Math. with Appl., 21.9 (1991): pp. 43-58.
% 

% written by Heather Wilber (heatherw3521@gmail.com)
% Jan, 2018. 

[m,n] = size(F); 
In = speye(n); 
Im = speye(m); 
compute_shifts = 0; 

%parse input 
if isempty(varargin)
    compute_shifts = 1;
    tol = eps; 
elseif numel(varargin)==1
    compute_shifts = 1; 
    tol = varargin{1}; 
elseif numel(varargin)==2
    p = varargin{1};
    q = varargin{2};
else %throw error
    error('ADI:adi: input invalid.')
end

% user wants shift parameters computed: 
if compute_shifts == 1
    %find intervals where eigenvalues live: 
    a = eigs(A, 1, 'SM'); 
    b = eigs(A, 1, 'LM'); 
    c = eigs(B, 1, 'SM'); 
    d = eigs(B, 1, 'LM'); 
    % determine if the eigenvalues have a complex part. if so, abandon. 
    if any(abs(imag([a b c d])) > 1e-10)
        error('ADI:adi:cannot automatically compute shift parameters unless the eigenvalues of A and B in AX - XB = F are contained in real, disjoint intervals.')
    end
    % check if intervals overlap
    I1 = [min(a,b) max(a,b)]; 
    I2 = [min(c,d) max(c,d)]; 
    if (( I1(1) < I2(2) && I1(2) > I2(1)) || ( I2(1) < I1(2) && I2(2) > I1(1)) )
        error('ADI:adi:cannot automatically compute shift parameters unless the eigenvalues of A and B in AX - XB = F are contained in real, disjoint intervals.')
    end
    [p,q] = getshifts_adi([I1 I2], 'tol', tol); 
end

%do ADI
X = zeros(m,n); 
for i = 1:length(p)
    X = (A - q(i)*Im)\(X*(B-q(i)*In)+F); 
    X = ( (A - p(i)*Im)*X - F)/(B-p(i)*In); 
end
