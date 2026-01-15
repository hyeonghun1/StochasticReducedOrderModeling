% function [E_test, C_test, f1, f2] = compute_model(f, V, x0_test, u_test, L_test)
function [E_test, C_test, f1, f2] = compute_model(f, V, xr0_test, s, L_test)

%Inputs--------------------------
% f: EulerMaruyamaStep function
% f = @(x0, u, L) EulerMaruyamaStep(x0, u, L, E, A, B, N, M, h)
% V: Vr POD basis (the identity R^{n x n} is used for the FOM).
% x0_test: I.C. of the testing data R^n
% u_test: input of the testing data R^s
% L_test: number of noise sampling, e.g., 1e6

%Outputs--------------------------
% E_test: mean of the FOM/ROM states R^{n x s} or R^{r x s}
% C_test: covariance of the FOM/ROM states (R^{n x n x s}) or R^{r x r x s}
% f1: second moment of X(end): ||X(end)||_2^2
% f2: mean over n final values X(T)^3 * exp(X(T))

r = size(xr0_test, 1);
% [m, s] = size(u_test);

% test parameters
batchSize = 170;  % Number of samples per batch
numBatches = ceil(L_test / batchSize);  % Total number of batches
% -> full simulation of 1,000,000 trajectories is performed with 100
% batches of 10,000 trajecs.

if L_test <= batchSize
  numBatches = 1;
  batchSize = L_test;
end

% initialize
E_test = zeros(r, s);
C_test = zeros(r, r, s);
f1 = 0;
f2 = 0;

% batch loop
for batch = 1:numBatches
  % batch
  % tic

  Nb = batchSize;
  if batch == numBatches
      Nb = L_test - (batch-1)*batchSize;  % last batch
  end

  % xr0_test = V'*xr0_test;

  % Compute stats for one batch
  % Each stat is sample resylt from one batch
  [E_temp, C_temp, f1_temp, f2_temp] = estimate(f, V, ...
                                            xr0_test, s, Nb);

  % Accumulate batch results
  E_test = E_test + (Nb/L_test)*E_temp;  % Works since all batches have the same size
  % C_test = C_test + 1/numBatches*C_temp;
  
  % Accumulate second moments E[X*X^T]
  for k=1:s
    C_test(:,:,k) = C_test(:,:,k) + (Nb/L_test) * ( C_temp(:,:,k) ...
                                                     + E_temp(:,k)*E_temp(:,k)' );
  end
  
  % f1 = f1 + 1/numBatches*f1_temp;
  % f2 = f2 + 1/numBatches*f2_temp;
  f1 = f1 + (Nb/L_test)*f1_temp;
  f2 = f2 + (Nb/L_test)*f2_temp;
  
  % toc; 
end

for k = 1:s
    C_test(:,:,k) = C_test(:,:,k) - E_test(:,k)*E_test(:,k)';
end

end