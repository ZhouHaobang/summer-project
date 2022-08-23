% Singular Value Shrinkage Operator

% Input
% T: threshold 
% A: a matrix waiting for shrinking singular values

% Output:
% X: the matrix after shringking singular values

function X = SVS(T,A)

% Singular value decomposition
[ U,S,V ] = svd(A);

% Thresholding
S = S - T;  

% Projection
SS = S > 0;
S = S.*SS;

% Recover the matrix 
X = U*S*V';
end
