% Test ADMM on a real image: Lena

% Read mage
I = imread('/Users/zhouhaobang/Desktop/Dissertation/Lenna.tif');

% Convert image to grayscale
I = rgb2gray(I);
I = im2double(I);

% size and rank of image
[m_h,n_h]=size(I);
rank_I = rank(I);

% Convert to low-rank image
[U,S,V] = svdsketch(double(I),1e-1);
Inew2 = U*S*V';
rank_Inew2 = rank(Inew2);
% imshow(Inew2)
% title(sprintf('Rank %d approximation',size(S,1)))

% Reconstruct on the original image
% B = I;

% Reconstruct on the low-rank image
B = Inew2;

% Hide some entries randomly
B_new = B;
for i=1:m_h
    for j=1:n_h
        p = rand;
        if p > 0.3
            B_new(i,j) = 0;
        end
    end
end

% Initialization
C = eye(m_h+n_h);
tol = 1e-4;
rho = 20;

% Projection
P = B_new > 0;
SR = sum(sum(P))/(m_h*n_h);  % Sampling Ratio

P_new = P';
P_new = P_new(:);

% Implement ADMM to reconstruct the image
tic
[X, W, primal, dual, gap, iter] = ADMM_SDP(C, B_new, P, P_new, m_h, n_h, tol, rho);
toc

% The relative error of the reconstructed image
rel_error = norm(W-B,'fro')/max(1,norm(B,'fro'));

% Show the reconstructed image
imshow(W)

