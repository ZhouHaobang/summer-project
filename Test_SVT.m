% Test the SVT algorithm on low-rank matrices

% Generate a low-rank matrix
[B, B_new, rank_B] = Create_low_rank_matrix(2, 100, 0.3);

% Parameters
[n1,n2] = size(B); 
T = 5*sqrt(n1*n2); % threshold to singular values
delta_t = 2.2;     % step size

% Projection
P = B_new > 0;

% SVT
tic
[ X,iterations,res,zhi] = SVT(B_new,P,T,delta_t);
toc

% Relative error
rel_error = norm(X-B,'fro')/max(1,norm(B,'fro'));

% X: the recovered matrix

