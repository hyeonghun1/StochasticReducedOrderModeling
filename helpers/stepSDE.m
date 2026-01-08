function X = stepSDE(f, x0, u, L)

%Inputs--------------
% f: EulerMaruyamaStep function
% f = @(x0, u, L) EulerMaruyamaStep(x0, u, L, E, A, B, N, M, h)
% x0: initial condition R^n
% u: input R^s
% L: number of noise samples

%Outputs--------------
% X: FOM observations in time, across noise realizations
% \in R^{n x L x s}
%---------------------

n = size(x0, 1);      % grid size
s = size(u, 2);       % time step
X = zeros(n, L, s);   % initialize snapshots

% First step
X(:,:,1) = f(x0, u(:,1), L);    % L samples of first snapshots (n x L)

% populate snapshots from second time step
for ii = 2:s
    X(:,:,ii) = f(X(:,:,ii-1), u(ii), L);
end

end