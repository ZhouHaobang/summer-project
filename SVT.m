% Singular value thresholding algorithm

% Input:
% A: an unknown matrix
% P: Projection
% T: a threshold
% delta_t: step size

% Output:
% X: the recovered matrix
% iter: the number of iterations
% res: residual error
% rank_X: rank of X


function [ X, iter, res, rank_X] = SVT(A, P, T, delta_t)

% Initialization
epsilon = Inf;
Y = zeros(size(A,1)); % An intermediate matrix
iter = 0;

while epsilon > 1e-4
       
        % Shrinks the singular values of matrix Y
        X = SVS(T,Y);

        % Update Y by the projection
        Y = Y + delta_t*P.*(A-X);
        
        % Compute the residual error
        epsilon = norm(P.*(A-X),'fro' )/norm(P.*A,'fro' )

        res(iter+1)=epsilon;
        rank_X(iter+1)=rank(X);
        
        iter = iter + 1
end
end


