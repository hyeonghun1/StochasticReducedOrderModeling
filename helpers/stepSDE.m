function Xr = stepSDE(f, x0, s, L)

%Inputs--------------
% f: EulerMaruyamaStep function
% f = @(x0, u, L) EulerMaruyamaStep(x0, u, L, E, A, B, N, M, h)
% x0: initial condition R^n
% s: total time steps
% L: number of noise samples

%Outputs--------------
% X: FOM observations in time, across noise realizations
% \in R^{r x L x s}
%---------------------

r = size(x0, 1);      % grid size
Xr = zeros(r, L, s);   % initialize snapshots

% First step
Xr(:,:,1) = f(x0, L);    % L samples of first snapshots (r x L)

% populate snapshots from second time step
for ii = 2:s
    Xr(:,:,ii) = f(Xr(:,:,ii-1), L);
end

end