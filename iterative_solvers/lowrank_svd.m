function [U, S, V] = lowrank_svd(C,R, varargin)
% computes the truncated SVD of C*R'.
%
% If C is an m by r matrix and R is an n by r matrix, then 
%
% lowrank_svd(C, R) computes the truncated SVD with U of size m by r, S of 
% size r by r, and R of size n by r. 
%
% lowrank_svd(C, R, p), where  0< p <= r, computes a rank p
% approximation to the SVD, so U and V each have p columns and S is p by p.
%
% lowrank_svd(C, R, 'tol', Tol) computes the truncated SVD of C*R', so that
%  min(diag(S)) > Tol*S(1,1) and length(S) <= r. 
%
 
% written by Heather Wilber (heatherw3521@gmail.com)
% Jan. 2018
%%
% compute SVD in low rank form: 

[Q1, R1] = qr(C,0); 
[Q2, R2] = qr(R,0); 

% Q1*R1*R2'*Q2' = C*R': 

[U, S, V] = svd(R1*R2'); 

if ~isempty(varargin)
    if numel(varargin)==3 %truncate to p sing vals.
        p = varargin{1}; 
        S = S(1:p, 1:p); 
        U = U(:,1:p); 
        V = V(:,1:p); 
    elseif nargin(varargin)==4 %truncate to tol
        p = varargin{2}; 
        idx = find(diag(S)>p*S(1,1), 1, 'last');
        S = S(1:idx, 1:idx); 
        U = U(:,1:idx); 
        V = V(:,1:idx); 
    end
end

U = Q1*U;
V = Q2*V; 

end


    

