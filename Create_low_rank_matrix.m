% Generate low-rank matrix

% Input:
% r: rank of the matrix
% x: size of the matrix
% sr: smapling ratio
% For example: Generate a 100-by-100 matrix of rank 2 with the 0.30 smapling ratio:
% [B, B_new, rank_B] = Create_low_rank_matrix(2,100, 0.3)

% Output:
% B: the original matrix
% B_new: undersampled matrix
% rank_B: rank of B

function [B, B_new, rank_B] = Create_low_rank_matrix(r, x, sr)

% Generate the original matrix
B = ones(1,x)' * ones(1,x) + rand(r-1,x)' * rand(r-1,x); 

% rank and size of B
rank_B=rank(B);      
[m_h,n_h]=size(B);

% Generate the undersampled matrix B_new
B_new = B;

% Hide entries according to the sampling ratio 'sr', and set them to zero
for i=1:m_h
    for j=1:n_h
        p = rand;
        if p > sr
            B_new(i,j) = 0;
        end
    end
end 

end
