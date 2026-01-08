% function [E_train, C_train, u_train] = train(f, Vr, m, s, L, h, isbilinear)

function [E_train, C_train, u_train] = train(f, Vr, m, s, L)

%Inputs -------------------------------
% f: EulerMaruyamaStep function
% f = @(x0, u, L) EulerMaruyamaStep(x0, u, L, E, A, B, N, M, h)
% Vr: POD basis of the total state snapshot X (n x L x s)
% m: input dimension
% s: number of time steps
% L: number of training noise realizations (e.g., (10, 100, 10000))

%Outputs ------------------------------
% E_train{j}{i}: mean of the reduced state for training case i
% using L{j} noise realizations
% C_train{j}{i}: covariance of the reduced state for training case i
% using L{j} noise realizations. Each element is R^{r x r}
% u_train{j}{i}: input signal used in training case i
%--------------------------------------

% The code generates training data of small condition number. This is done
% with two parts below (zero initial condition/constant inputs
%                       & zero inputs/I.Cs using dominant POD modes)


L_max = max(L);

[n, rmax] = size(Vr);

% Generate training data
numL = numel(L);
E_train = cell(numL, 1);    % empirical mean of reduced states     
C_train = cell(numL, 1);    % empirical covariance of reduced states
u_train = cell(numL, 1);    % input signals

% To generate a data-matrix of small condition number, we repeatedly
% evaluate the FOM over pais of I.C.s and input.
% Part 1 - Training with different constant inputs
% This is done by querying FOM using zero I.C.s and constant inputs
k = 21;
us = linspace(-2,2,k);  % constant input amplitudes

idx = 0;
for ii=1:k
    idx = idx+1;
    disp("x0 = 0, u = e_"+ii);
    x0 = zeros(n,1);        % zero initial state
    u = us(ii)*ones(m, s);  % constant-in-time control input

    % FOM trajectory with L_max stochastic realizations                        
    X_temp = stepSDE(f, x0, u, L_max);  % R^{n x L_max x s}
    
    % Project the FOM observations
    X_train = pagemtimes(Vr', X_temp);  % R^{r x L_max x s}

    % For each number of noise sampling,
    % estimate mean/covariance of the (reduced) projected observations.
    for jj = 1:numel(L)
        LL = L(jj);     % 10 or 100 or 10000
        
        % mean of reduced states across LL realizations
        E_train{jj}{idx} = squeeze(mean(X_train(:, 1:LL, :), 2));

        % covariance of reduced states
        C_train{jj}{idx} = page_cov(X_train(:, 1:LL, :), true) ;
        
        % Store input u
        u_train{jj}{idx} = u;
    end
end

% Part 2 - Additional reduced training data (each POD mode as initial
% conditions)
% This is done by querying FOM r times more.
for ii=1:rmax
    idx = idx + 1;
    disp("x0 = Vr(:,"+ii+"), u = 0");
    x0 = Vr(:,ii);     % initial conditions using dominant POD modes
    u = zeros(m,s);    % zero control inputs

    X_temp = stepSDE(f, x0, u, L_max);
    X_train = pagemtimes(Vr', X_temp);
    for jj = 1:numel(L)
        LL = L(jj);
        E_train{jj}{idx} = squeeze(mean(X_train(:, 1:LL, :),2));
        C_train{jj}{idx} = page_cov(X_train(:, 1:LL, :), true) ;
        u_train{jj}{idx} = u;
    end
end


end


