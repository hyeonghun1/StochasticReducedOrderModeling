function [E_emp, C_emp, f1, f2] = estimate(f, Vr, xr0, s, L)

%Inputs -------------------------------
% f: EulerMaruyamaStep function
% f = @(x0, u, L) EulerMaruyamaStep(x0, u, L, E, A, B, N, M, h)  
% Vr: POD basis of the total snapshot X (n x L x s)
% xr0: initial condition R^r
% s: total time steps
% L: number of noise samples

%Outputs -------------------------------
% E: mean of the reduced states R^{n x s}
% C: covariance of the reduced states Vr*C*Vr^T (R^{n x n x s})
% f1: second moment of X(end): ||X(end)||_2^2
% f2: mean over n final values X(T)^3 * exp(X(T))
%---------------------------------------

% Run L stochastic simulations starting from x0 with input u
Xr = stepSDE(f, xr0, s, L);   % R^{r x L x s}

% Estimate quantities
% Empirical mean over L noise realizations
E_emp = reshape(mean(Xr,2), size(Xr, [1 3])); % R^{r x s}

% E = Vr * E_emp;  % mean R^{n x s}  

% Empirical covariance 
C_emp = page_cov(Xr, true);

% C = pagemtimes( pagemtimes(Vr, C_emp), Vr' );  % Vr*C*Vr^T (R^{n x n x s})

% Compute diagonostics to measure wack errors at the final time.
X_T_recon = Vr * Xr(:,:,end);  % final state reconstruction for all L noises (R^{r x L})
f1 = mean(vecnorm(X_T_recon.^2));        
f2 = mean(X_T_recon.^3./exp(X_T_recon), 'all');

end

