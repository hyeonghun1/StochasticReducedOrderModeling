function [Mhat, Khat] = infer_diffusion(C_train, h, Ahat, Nhat)

%Inputs--------------------------
% C_train: each cell contains the covariance of the reduced states
% C_r(t) = cov(Xr, Xr) = Vr^T * cov(X, X) * Vr -> R^{r x r}
% u_train: input matrix 
% h: sampling time dt
% Ahat: linear reduced operator, inferred from drift OpInf
% Nhat: bilinear reduced operator, inferred from drift OpInf

%Outputs--------------------------
% Mhat: reduced diffustion matrix 
% Khat: correlation matrix (this is set to identity)
% -> Hhat = Mhat * Khat * Mhat^T (diffusion contribution)
%---------------------------------

% m = size(u_train{1}, 1);    
% p = numel(u_train);

% The solution of the least squares problem is the mean over the residuals 
Hhat = zeros(size(Ahat));

[Cr, Cr_dot, ind] = central_finite_differences(C_train, h, 2, 3);
% u = u_train{ii};

for jj=1:numel(ind)
    % Clyap: \Psi_r(t) * Cr(t) in the paper..
    Psi_hat = Ahat;
    Clyap = Psi_hat * Cr(:,:,jj);
    
    % Update Hhat (diffusion contribution)
    incre = Cr_dot(:,:,jj) - Clyap - Clyap';
    % Hhat = Hhat + 1/(p*numel(ind)) * incre; 
    Hhat = Hhat + 1/numel(ind) * incre;
end


% Make sure that the Hhat is symmetric
Hhat = Hhat/2 + Hhat'/2;

% Eigendecomposition of H and sort
[HU, HS] = eig(Hhat);
[Ssort, I] = sort(diag(HS), 'descend');
HS = diag(Ssort);
HU = HU(:, I);

% get the index of the last singular value that is greater or equal
% than tol*\sigma_max: (corresponds to 2-norm)
tol = max(HS, [], 'all')/1000;  % a threshold of 0.1% of the max eigenvalue
d = find(diag(HS) >= tol, 1, 'last');
if isempty(d)
  d=0;
end

% truncate H and decompose
Mhat = HU(:,1:d) * sqrt(HS(1:d,1:d));
Khat = eye(d);   % corrrelation matrix of the surrogate noise process
end

