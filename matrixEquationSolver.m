function X = matrixEquationSolver(A, B, F)
% MATRIXEQUATIONSOLVER 
%
% X = matrixEquationSolver(A, B, F) solves the generalized Sylvester
% equation: 
%
%      sum_j A{j} X B{j}^T = F, 
%
% where A and B are cell array of nxn and mxm matrices, respectively. 
% 
% There is no known efficient and general algorithm when the number of terms
% in the sum is more than 2. The resulting algorithm based on Kronecker 
% products costs O((mn)^3) operations.  

% Written by Alex Townsend, Nov 2014. (alex.townsend1987@gmail.com)

% Do not form matrices above this size: 
maxSize = 5000; 

if ( ~isa(A, 'cell') || ~isa(B, 'cell') ) 
    error('MATRIXEQUATIONSOLVER:INPUTS:CELL',...
                      'Defining matrices are not in a cell array.'); 
end

% Check that the sizes of the cell arrays are consistent: 
if ( size(A) ~= size(B) ) 
    error('MATRIXEQUATIONSOLVER:INPUTS:SQUARE', ...
                       'Ambigous number of terms in the matrix equation.'); 
end 

% All the matrices should be square and of the same size: 
[mA, nA] = size( A{ 1 } ); 
[mB, nB] = size( B{ 1 } ); 
if ( (mA ~= nA) || (mB ~= nB) )
    error('MATRIXEQUATIONSOLVER:INPUTS:SIZES', ...
                           'Rectangular matrices are not allowed.'); 
end

% If we have two terms, then use bartelsStewart: 
if  ( ( max(size(A)) == 1 ) )
 % Solving A X B^T = F
    Y = A{1} \ F;    % Y = XB^T
    X = Y / B{1}.';  % Solve for X
elseif ( ( max(size(A)) == 2 ) )
    X = bartelsStewart(A{1}, B{1}, A{2}, B{2}, F ); 
else
    % Form the large Kronecker matrix. 
    
    % Check that the matrices are not too large: 
    if ( nA * nB > maxSize ) 
        error('MATRIXEQUATIONSOLVER:MEMORY', ... 
                                  'A very large matrix will be formed.' );
    end

    % Form the large nA x nB matrix: 
    C = zeros( nA*nB ); 
    for j = 1:max(size(A))
        C = C + kron( B{j}, A{j} ); 
    end
    
    % Solve the large linear system: 
    x = C \ F(:); 
    
    % Reshape solution vector to a matrix: 
    X = reshape( x, nA, nB ); 
end

end