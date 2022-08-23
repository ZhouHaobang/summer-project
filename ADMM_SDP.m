% Implement ADMM for the SDP reformulation of matrix completion

% Dual problem
% min  b^Ty
% s.t. Ay + S = C
%      S >= 0

% SDP reformulation
% min  <I_n, X>
% s.t. <A_L, X> = b_l, l=1,...,m
%      X >= 0

% Input:
% C: A given matrix
% B_new: an undersampled matrix
% P: a projection matrix such that P(i,j)=1 if B_new(i,j)!=0, P(i,j)=0, otherwise.
% P_new: the column vector of P
% m_h: the number of rows of the matrix B_new
% n_h: the number of columns of the matrix B_new
% tol: tolerance
% rho: a penalty parameter

% Output: 
% X: the variable of the SDP
% W: the recovered matrix
% primal, dual, gap: three different infeasibilities
% iter: the number of iterations

function [X, W, primal, dual, gap, iter] = ADMM_SDP(C, B_new, P, P_new, m_h, n_h, tol, rho)

% Size of matrix X
m = m_h*n_h;
n = m_h+n_h;

% Transform the know matrix B to a column vector b
B_ = B_new';
b = B_(:);

% Initialization
delta = Inf;
X = eye(n);
S = eye(n);
iter = 0;

% Parameters initialization for updating the penalty parameter
gamma = 0.5;
counter_d = 0;
counter_i = 0;

% Initially compute AX
AX = zeros(m,1);
k = 0;
for i = 1:m_h
     for j = 1:n_h
         k = k + 1;
         if P(i,j) == 1
             AX(k) = X(i,m_h+j) + X(m_h+j,i);
         else
             AX(k) = 0;
         end
     end    
end

% Iterations
while (delta > tol)

    iter = iter + 1 % Record the number of iterations
   
    % Compute A(S-C)
    t = 0;
    A_SC = zeros(m,1);
    for i = 1:m_h
        for j = 1:n_h
            t = t + 1;
            if P(i,j) == 1
                A_SC(t) = S(i,m_h+j) + S(m_h+j,i) - C(i,m_h+j) - C(m_h+j,i);
            else
                A_SC(t) = 0;
            end
        end    
    end
        
    % Update y
    y = -(0.5).*((1/rho).*(AX-b) + A_SC);
    
    % Compute A_star_y
    W_update = P_new.*y;
    W_update = reshape(W_update,m_h,n_h)';
    A_star_y = [zeros(m_h,m_h),W_update; W_update',zeros(n_h,n_h)];
    
    
    % Update Z
    Z = C - A_star_y - (1/rho).* X;
 
    % eigenvalues decomposition, D is a diagonal matrix with eigenvalues
    [V, D] = eig(Z);
    [L,LD] = ldl(Z);
    con(iter) = cond(L*L');
    
    % Store the eigenvalues before the projection 
    Eig_B(:,iter) = diag(D);
    
    % Projection: cut the negative eigenvalues
    DD = D >= 0;
    D = DD.*D;
    
    % Store the eigenvalues after the projection
    Eig_A(:,iter) = diag(D);
    
    % Update S
    S = V*D*V';
    
    % Update X
    X = rho.*(S - Z);
    
    % Update AX
    AX = zeros(m,1);
    k = 0;
    for i = 1:m_h
        for j = 1:n_h
            k = k + 1;
            if P(i,j) == 1
                AX(k) = X(i,m_h+j) + X(m_h+j,i);
            else
                AX(k) = 0;
            end
        end    
    end
    
    % Primal infeasibility
    primal_up = norm(AX - b, 2);
    primal_down = norm(b,2);
    primal = primal_up/primal_down;
%     primal = norm(AX - b, 2)/ norm(b,2);
    pinf(iter) = primal; % Store the primal infeasibility

    % Dual infeasibility
    dual_up = (norm(A_star_y + S - C, 'fro'))^2;
    dual_down = (1+norm(C,1));
    dual = dual_up/dual_down;
%     dual = (norm(A_star_y + S - C, 'fro'))^2/(1+norm(C,1));
    dinf(iter) = dual;  % Store the dual infeasibility

    % Gap infeasibility
    gap_up = abs(b'*y - trace(C'*X));
    gap_down = (1+abs(b'*y)+trace(C'*X));
    gap = gap_up/gap_down;
%     gap = abs(b'*y - trace(C'*X))/(1+abs(b'*y)+trace(C'*X));
    ginf(iter) = gap;
    
    delta = primal
%     delta = max([primal, gap, dual])
    
    
    % Updating the Penalty Parameter if m_h > 100
    if primal_up < dual
        counter_d = counter_d + 1;
    else
        counter_d = 0;
    end
    
    if counter_d > 10
        increment = rho/gamma;
        if increment > 1e-4 && increment < 1e+4
            rho = increment;
        end    
    end
    
    if primal_up > dual
        counter_i = counter_i + 1;
    else
        counter_i = 0;
    end
    
    if counter_i > 10
        decreasing = rho*gamma;
        if decreasing > 1e-4 && decreasing < 1e+4
            rho = decreasing;
        end    
    end    
    
    % Recovererd matrix W
    W = 2*X(1:m_h,m_h+1:m_h+n_h);
%     r_W(iter) = rank(W);
    
end
end
