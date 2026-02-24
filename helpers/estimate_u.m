function [E_emp, C_emp] = estimate_u(f, Vr, xr0, u, L)

%Inputs -------------------------------
% f: EulerMaruyamaStep function
% f = @(x0, u, L) EulerMaruyamaStep(x0, u, L, E, A, B, N, M, h)  
% Vr: POD basis of the total snapshot X (n x L x s)
% xr0: initial condition R^r
% u: input
% L: number of noise samples

%Outputs -------------------------------
% E: mean of the reduced states R^{n x s}
% C: covariance of the reduced states Vr*C*Vr^T (R^{n x n x s})
% f1: second moment of X(end): ||X(end)||_2^2
% f2: mean over n final values X(T)^3 * exp(X(T))
%---------------------------------------

% Run L stochastic simulations starting from x0 with input u
Xr = stepSDE_u(f, xr0, u, L);   % R^{r x L x s}

% Estimate quantities
% Empirical mean over L noise realizations
E_emp = squeeze(mean(Xr, 2));  % R^{r x s}
% E_emp = reshape(E_emp, size(Xr, [1 3])); % R^{r x s}

% Empirical covariance 
C_emp = page_cov(Xr, true);

% E = Vr * E_emp;  % mean R^{n x s}
% C = pagemtimes( pagemtimes(Vr, C_emp), Vr' );  % Vr*C*Vr^T (R^{n x n x s})

end

