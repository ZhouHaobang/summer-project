% A modified Cholesky algorithm

% P(A+E)P^T = L*DMC*L^T  

% Input:
% A: a symmetric indefinite matrix
% delta: a threshold. If we want a PSD matrix, we set delta = 0.

% Output:
% L: a unit lower triangular matrix
% DMC: a block diagonal matrix with diagonal blocks of 1 or 2 after the modification
% P: a permutation matrix
% D: a block diagonal matrix with diagonal blocks of 1 or 2
% A_new: a positive definite (positive semidefinite) matrix

% Example: a symmetric indefinite matrix for testing
% A = [6 12 3 -6; 12 -8 -13 4; 3 -13 -7 1; -6 4 1 6];

function [L, DMC, P, D, A_new] = modchol_ldlt(A, delta)

% Determine whether the matrix is symmetric
if ~ishermitian(A)
    error('Must supply symmetric matrix.')
end

% if we did not give the value of delta, then the dalta is
if nargin < 2
    delta = sqrt(eps)*norm(A,'fro'); 
end

n = max(size(A));

% Cholesky factorization for an indefinite matrix
[L,D,p] = ldl(A,'vector'); 

% Initialization
DMC = eye(n);

% Modified Cholesky perturbations.
k = 1;
while k <= n

      if k == n || D(k,k+1) == 0 % 1-by-1 block
         
         % Threshold delta 
         if D(k,k) <= delta
            DMC(k,k) = delta;
         else
            DMC(k,k) = D(k,k);
         end
         k = k+1;
      
      else % 2-by-2 block
         
         % Threshold delta 
         E = D(k:k+1,k:k+1);
         [U,T] = eig(E);
         for ii = 1:2
             if T(ii,ii) <= delta
                T(ii,ii) = delta;
             else
                T(ii,ii) = T(ii,ii); 
             end
         end
         
         temp = U*T*U'; % new PSD matrix E
         DMC(k:k+1,k:k+1) = (temp + temp')/2;  % Guarantee symmetric.
         k = k + 2;

      end

end

% Permutation matrix
P = eye(n); 
P = P(p,:); 

% The positive definite (semidefinite) matrix 
A_new = P*L*DMC*L'*P';

end

%   Reference:
%   S. H. Cheng and N. J. Higham. A modified Cholesky algorithm based
%   on a symmetric indefinite factorization. SIAM J. Matrix Anal. Appl.,
%   19(4):1097-1110, 1998. doi:10.1137/S0895479896302898,

%   Authors: Bobby Cheng and Nick Higham, 1996; revised 2015.
%   Revised by Haobang Zhou