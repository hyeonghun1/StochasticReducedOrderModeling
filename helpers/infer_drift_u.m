function [Ehat,Ahat,Bhat,Nhat] = infer_drift_u(E_train, u_train, h, isbilinear, lambda)

%Inputs--------------------------
% E_train: Expected value of the reduced states across samples R^{r x s}
% E_r(t) = E[Vr^T X(t)]
% u_train: input matrix 
% h: dt
% lambda: pernalizing term

%Outputs--------------------------
% Ehat: mass matrix
% Ahat: Ar_hat (linear reduced operator)
% Bhat: Br_hat (input reduced operator)
% Nhat: Nr_hat (bilinear reduced operator)
%---------------------------------

r = size(E_train, 1);

[m, s] = size(u_train);
u = u_train;

% central finite difference quotient of accuracy 2 along dimension 2
[Er, Er_dot, ind] = central_finite_differences(E_train, h, 2, 2);

if isbilinear
    col_dim = r + m + r*m;
else
    col_dim = r + m;
end

K = min(numel(ind), s);

D   = zeros(col_dim, K);
rhs = zeros(r, K);

if isbilinear
    % Needs kronecker product btw. input and Er
    for jj = 1:K
        D(:,jj) = [Er(:,jj); u(:,jj); kron(u(:,jj), Er(:,jj))];
        rhs(:,jj) = Er_dot(:,jj);
    end
else
    D = [Er(:,1:K); u(:,1:K)];
    rhs = Er_dot(:,1:K);
end

% fprintf("cond(D) = %.4e\n", cond(D))

% Solve least-squares problem
if exist('lambda', 'var')
    % Get regularization matrix
    Gamma = regularizer(r, lambda, isbilinear);
    
    % Solve regulairzed least-squares problem
    Dt = D' ; Rt = rhs';
    D_modified = Dt' * Dt + Gamma'*Gamma;
    rhs_modified = D * Rt;
    ops_t = D_modified \ rhs_modified;
    ops = ops_t';
    Ahat = ops(1:r, 1:r);
    Bhat = ops(:, r+1:r+m);
else
    % Non-regularized least-squares problem
    Dt = D'; Rt = rhs';
    ops = (Dt \ Rt)';
    % ops = rhs/D;
    Ahat = ops(1:r, 1:r);
    Bhat = ops(:, r+1:r+m);
end

if isbilinear
    Nhat = ops(:, r+m+1:end);
else 
    Nhat = zeros(size(Ahat));
end

% In our examples E is always the identity
Ehat = speye(r);
end

