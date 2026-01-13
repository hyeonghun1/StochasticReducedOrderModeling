function [Ehat,Ahat,Nhat] = infer_drift(E_train, h, isbilinear, s, reg)

%Inputs--------------------------
% E_train: Expected value of the reduced states across samples R^{r x s}
% E_r(t) = E[Vr^T X(t)]
% u_train: input matrix 
% h: dt
% doReg: logical (true/false) for regularization

%Outputs--------------------------
% Ehat: mass matrix
% Ahat: Ar_hat (linear reduced operator)
% Bhat: Br_hat (input reduced operator)
% Nhat: Nr_hat (bilinear reduced operator)
%---------------------------------

r = size(E_train, 1);
D = [];
rhs = [];

% m = size(u_train{1}, 1);
% p = numel(E_train);

% central finite difference quotient of accuracy 2 along dimension 2
[Er, Er_dot, ind] = central_finite_differences(E_train, h, 2, 2);

% u = u_train{ii};

for jj=1:min(numel(ind), s)
    % disp(jj)

    if isbilinear
      % system has bilinear term
      D_new = [D, [Er(:,jj); u(:,jj) kron(u(:,jj), Er(:,jj))] ];
    else
      % system does not have a bilinear term
      D_new = [D, [Er(:,jj)]];
    end

    % if ii<2 %r+m
        % use full trajectories
        D = D_new;
    % else
        % only use this time point if it reduces the condition number
    %     if cond(D_new) < cond(D)
    %         D = D_new;
    %     else
    %         continue
    %     end
    % end

    rhs = [rhs, Er_dot(:,jj)];
end

disp("cond(D) = " + cond(D))


if exist('lambda', 'var') && ~isempty(reg)
    % Get regularization matrix
    Gamma = regularizer(r, reg);
    
    % Solve regulairzed least-squares problem
    D_modified = D'*D + Gamma'*Gamma;
    rhs_modified = D' * R';
    ops = rhs_modified/D_modified;
    Ahat = ops(1:r, 1:r);
else
    % Non-regularized least-squares problem
    ops = rhs/D;
    Ahat = ops(1:r, 1:r);
    % Bhat = ops(:, r+1:r+m);
end

if isbilinear
    Nhat = ops(:, r+m+1:end);
else 
    Nhat = zeros(size(Ahat));
end

% In our examples E is always the identity
Ehat = speye(r);
end

