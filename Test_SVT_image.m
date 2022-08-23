% Test the SVT algorithm on a real image

% Read image
I = imread('/Users/zhouhaobang/Desktop/Dissertation/Lenna.tif');

% Convert image to grayscale
I = rgb2gray(I);

% Convert to low-rank image
[U,S,V] = svdsketch(double(I),1e-1);
Inew2 = U*S*V';
rank_Inew2 = rank(Inew2);

% Size of the image
[m_h,n_h]=size(Inew2);

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

% Parameters
[n1,n2] = size(B); 
T = 5*sqrt(n1*n2); % threshold to singular values
delta_t = 1.2;     % step size

% Projection
P = B_new > 0;

% SVT
tic
[ X,iter,res,rank_X] = SVT(B_new,P,T,delta_t);
toc

% Relative error
rel_error = norm(X-B,'fro')/max(1,norm(B,'fro'));

% Uint8 convert
X = uint8(X);

% Show the reconstructed image
imshow(X)