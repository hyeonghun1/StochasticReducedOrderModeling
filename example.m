clear all; clc;
addpath(genpath("./"))
rng(2)

oneDflag = true;

if oneDflag
    [E, A, B, N, M, f, n, m, h, s] = heatExample();
else
    [E, A, B, N, M, f, n, m, h, s, ind] = heat2dExample();
end

% Bilinearity check
isbilinear = full(any(N, 'all'));

rmax = 10;  % maximum ROM dimension

% number of samples
L_subspace = 1e4;
L_train = [1e1 1e2 1e4];
L_test = 1e6;

%% Plot one trajectory of 1D/2D Heat SDE
if oneDflag
    % 1D Heat SDE solution
    dx = 1/(n+1);
    x = (1:n)*dx;      % spatial grid
    Nt = 1000;         % Number of time steps
    x0 = zeros(n,1);   % Zero initial condition
    
    % % V-shaped spline connected to Dirichlet BC
    % u = 1;  % Dirichlet BC value at boundaries
    % x0 = u - 2*u*min(x, 1-x);
    
    % Dirichlet BC: u(t) = cos(2*pi*t/T_total)
    u_subspace = cos(linspace(0,2*pi,Nt));
           
    % Initialize storage
    X = zeros(n, Nt);
    X(:,1) = x0;
    
    % SDE time-stepping
    for k = 2:Nt
        X(:,k) = f(X(:,k-1), u_subspace(k), 1);  % advance one step
    end

    % figure;
    % plot(x, X(:,3), 'b', 'LineWidth', 2);
    
    % Animate the 1D heat solution over time
    figure;
    hLine1 = plot(x, X(:,1), 'b', 'LineWidth', 2);
    xlabel('x'); ylabel('y(x,t)');
    ylim([min(X(:)) max(X(:))]);
    grid on;
    title(sprintf('1D Heat SDE: t = %.3f', 0));
    
    % Loop over time steps
    for k = 1:Nt
        set(hLine1, 'YData', X(:,k));
        title(sprintf('1D Heat SDE: t = %.3f', h*k));
        drawnow;
    end

else
    % 2D Heat SDE solution
    nx = 34;
    Nt = 100;
    X_full = zeros(nx, nx, Nt);
    
    % Dirichlet BC: u(t) = cos(2*pi*t/T_total)
    u_subspace = cos(linspace(0,2*pi,Nt));
    x0_subspace = zeros(n, 1);
    
    % Initialize storage
    X_active = zeros(n, Nt);
    
    % SDE time-stepping
    for k = 2:Nt
        X_active(:,k) = f(X_active(:,k-1), u_subspace(k), 1);  % advance one step
    end

    % Reshape each vector (active) to 2D full grid
    for k = 1:Nt
        tmp = zeros(nx^2,1);
        tmp(ind) = X_active(:,k);            
        X_full(:,:,k) = reshape(tmp, nx, nx);
    end
    
    % 2D animation
    fig = figure;
    hImg = imagesc(X_full(:,:,1));  % initial frame
    axis equal tight;
    colorbar;
    caxis([min(X_active(:)) max(X_active(:))]);
    xlabel('x'); ylabel('y');
    title(sprintf('2D Heat SDE'));
    
    frame_idx = 1;
    fps = 10;            % desired frames per second
    pause_time = 1/fps;  % time between frames
    
    while ishandle(fig)        % loop until figure is closed
        set(hImg, 'CData', X_full(:,:,frame_idx));
        title(sprintf('2D Heat SDE, Time step %d', frame_idx));
        drawnow;
   
        frame_idx = frame_idx + 1;  % advance frame
        if frame_idx > Nt
            frame_idx = 1;          % repeat from start
        end
        pause(pause_time);          % slow down
    end

    % 3D Surface animation and export to MP4

% 3D Surface animation with fixed surf z-axis

fig = figure;

[nx, ny, Nt] = size(X_full);  % dimensions
[xgrid, ygrid] = meshgrid(1:nx, 1:ny);

% Compute global z-limits for surf
zmin = min(X_active(:));
zmax = max(X_active(:));

% Initial surface
hSurf = surf(xgrid, ygrid, X_full(:,:,1));
shading interp
colormap turbo
colorbar
caxis([zmin zmax]);       % fixed color scale
zlim([zmin zmax]);        % fixed z-axis
xlabel('x'); ylabel('y'); zlabel('u(x,y)');
title('2D Heat SDE (3D Surface)');
view(45, 35);
axis tight

fps = 10;

% Create VideoWriter object
v = VideoWriter('Heat2D_Animation.mp4','MPEG-4');
v.FrameRate = fps;
open(v);

% Loop through frames
for frame_idx = 1:Nt
    set(hSurf, 'ZData', X_full(:,:,frame_idx));  % update surface height
    set(hSurf, 'CData', X_full(:,:,frame_idx));  % update surface color
    title(sprintf('2D Heat SDE, Time step %d', frame_idx));
    caxis([zmin zmax]);       % fixed color scale
    zlim([zmin zmax]);        % fixed z-axis
    drawnow;
    colormap jet;
    % Capture the frame
    frame = getframe(fig);
    writeVideo(v, frame);
end

close(v);

disp('MP4 video saved as Heat2D_Animation.mp4');



end


%% Collect observations
u_subspace = cos(linspace(0,2*pi,s));  % time-dependent BC (input)
x0_subspace = zeros(n,1);       % zero initial condition

disp("Collecting all FOM observations X(t_i, w_j), i=1,...,s, j=1,...,L")
X_subspace = stepSDE(f, x0_subspace, u_subspace, L_subspace);
disp("State snapshot computation completed.")


%% Animate one trajectory of the simulated 1D Heat SDE solution

sample_idx = 2;  % first stochastic sample

X1 = X_subspace(:,sample_idx,:);  % space x time

[nx, Nt] = size(X1);
x = linspace(0, 1, nx);          % spatial grid
t = linspace(0, (Nt-1)*h, Nt);   % time vector

figure;
hLine1 = plot(x, X1(:,1), 'LineWidth', 1.5);
xlabel('Space');
ylabel('Amplitude');
ylim([min(X1(:)) max(X1(:))]);
grid on;
title(sprintf('1D Heat equation: Sample %d, t = %.3f', ...
        sample_idx, t(1)));

for k = 1:Nt
    set(hLine1, 'YData', X1(:,k));
    title(sprintf('1D Heat equation: Sample %d, t = %.3f', ...
            sample_idx, t(k)));
    drawnow;
end


%% Animate multiple samples of 1D Heat SDE solution

sample_idx = [1, 2, 3, 10, 20, 30, 100, 200, 300];  % indices of samples to animate
num_samples = length(sample_idx);

[nx, ~, Nt] = size(X_subspace);  % X_subspace: space x sample x time
x = linspace(0, 1, nx);          % spatial grid
t = linspace(0, (Nt-1)*h, Nt);   % time vector

figure;
hold on;
hLine1 = gobjects(1, num_samples);

colors = lines(num_samples);       % choose distinct colors
for ts = 1:num_samples
    Xs = squeeze(X_subspace(:, sample_idx(ts), :)); % nx x Nt
    hLine1(ts) = plot(x, Xs(:,1), 'LineWidth', 1.5, 'Color', colors(ts,:));
end

xlabel('Space');
ylabel('Amplitude');
ylim([min(X_subspace(:)) max(X_subspace(:))]);
grid on;

% title('1D Heat equation: Multiple samples');

% Animation loop
while true
    for k = 1:Nt
        for ts = 1:num_samples
            Xs = squeeze(X_subspace(:, sample_idx(ts), :)); % nx x Nt
            set(hLine1(ts), 'YData', Xs(:,k));
        end
        title(sprintf('1D Heat equation: t = %.3f', t(k)));
        drawnow;
    end
end


%% Animate one trajectory of the simulated 2D Heat SDE

sample_idx = 100;
X1 = X_subspace(:, sample_idx, :);   % n_active x Nt
[n_active, Nt] = size(X1);
n = 34;
X_full = zeros(n, n, Nt);

% Reshape each vector (active) to 2D full grid
for k = 1:Nt
    tmp = zeros(n^2,1);
    tmp(ind) = X1(:,k);            
    X_full(:,:,k) = reshape(tmp, n, n);
end

% Continuous 2D Animation (slower)
fig = figure;
hImg = imagesc(X_full(:,:,1));    % initial frame
axis equal tight;
colorbar;
caxis([min(X1(:)) max(X1(:))]);    % fixed color scale
xlabel('x'); ylabel('y');
title(sprintf('2D Heat Equation Sample %d', sample_idx));

frame_idx = 1;
fps = 10;                 % desired frames per second
pause_time = 1/fps;       % time between frames

while ishandle(fig)        % loop until figure is closed
    set(hImg, 'CData', X_full(:,:,frame_idx));
    title(sprintf('2D Heat Equation Sample %d, Time step %d', sample_idx, frame_idx));
    drawnow;

    % advance frame
    frame_idx = frame_idx + 1;
    if frame_idx > Nt
        frame_idx = 1;     % repeat from start
    end

    pause(pause_time);     % slow down
end


%% Get POD basis

if ~exist("V", "var")  
    disp("Computing subspace Vr...")
    [V,S,~] = svd(reshape(X_subspace, n, []), "econ");
end
disp("Subspace computation done!")

Vrmax = V(:, 1:rmax);

%%
normalized_sigma = diag(S)/max(diag(S));
figure;
semilogy(normalized_sigma);
hold on;
scatter(1:length(diag(S)), normalized_sigma, 'filled');hold off;
xlabel('Index');
ylabel('Singular Values');
grid on;


%% Get the ROM training data
% This generates mean/covariance of the reduced trajectories across
% all noise realizations, which are E_train, C_train.
% In the paper, this is denoted as E_r^L and C_r^L, which are needed for
% the data in drift and diffusion OpInf. The input u_train is also needed
% in the data matrices.

[E_train, C_train, u_train] = train(f, Vrmax, m, s, L_train);

% E_train/C_train has three elements (3 cases of noise), each of which
% has (k+r) pairs. The FOM is repeatedly queried to obtain a data-matrix
% of small condition number.
% Each (10 x 1000)/(10 x 10 x 1000) is the trajectory of 
% the expected/covariance values of the reduced states.


%% test set
u_test = cos(linspace(0, 5*pi, s)); % control input for testing
x0_test = x0_subspace;    % zero initial condition
seed_test = 42;

% FOM reference
disp("Computing FOM reference (testing) data...")
rng(seed_test)

% For the FOM (testing), identity matrix is used for Vr
[EFOM, CFOM, f1FOM, f2FOM] = compute_model( ...
    f, eye(n), x0_test, u_test, L_test);

disp("FOM tesing computed!")


%% plot FOM Expectation (EFOM)

% FOM expectation
[nx, Nt] = size(EFOM);
x = linspace(0, 1, nx);          % spatial grid
t = linspace(0, (Nt-1)*h, Nt);   % time vector

figure()
hLine1 = plot(x, EFOM(:,1), 'LineWidth', 1.5);
xlabel('Space');
ylabel('Amplitude');
ylim([min(EFOM(:)) max(EFOM(:))]);
grid on;
title(sprintf('FOM Expectation 1D Heat equation t = %.3f', t(1)));

for k = 1:Nt
    set(hLine1, 'YData', EFOM(:,k));
    title(sprintf('FOM Expectation 1D Heat equation t = %.3f', t(k)));
    drawnow;
end

%% plot FOM covariance
Nt = size(CFOM, 3);
fps = 50;                 % desired frames per second
pause_time = 1/fps;       % time between frames

figure;
im = imagesc(CFOM(:,:,1))
colormap(jet);
colorbar;
caxis([min(CFOM(:)) max(CFOM(:))]);

for kk=1:Nt
    set(im, 'CData', CFOM(:,:,kk));
    title(sprintf('FOM Covariance 1D Heat t = %.3f', t(kk)));
    drawnow;
    pause(pause_time);
end


%% OpInf
% We infer reduced operators, simulate ROM, get moments, and
% compute errors. We also test three different noise samplings.
E_error_OpInf = zeros(rmax, numel(u_train));
C_error_OpInf = zeros(rmax, numel(u_train));
f1_error_OpInf = zeros(rmax, numel(u_train));
f2_error_OpInf = zeros(rmax, numel(u_train));

% Infer using different number of noise samples
for kk = 1:numel(E_train)
    % disp(kk)
    [Ehat, Ahat, Bhat, Nhat] = infer_drift(E_train{kk}, u_train{kk}, h, isbilinear, s);
    [Mhat, Khat] = infer_diffusion(C_train{kk}, u_train{kk}, h, Ahat, Nhat);

    % Vary ROM dimension from 1 to rmax (for plotting)
    for ii=1:rmax
      Vr = Vrmax(:, 1:ii);
      Ehatr = Ehat(1:ii, 1:ii);
      Ahatr = Ahat(1:ii, 1:ii);
      Bhatr = Bhat(1:ii, :);
      Nhatr = Nhat(1:ii, 1:ii);
      Mhatr = Mhat(1:ii, :);

      % ROM propagation via EM time step, using inferred operators
      % x_{k+1} = (E_r - h*A_r - h*N_r*u)^-1 * [x_k + h*B_r*u + \sqrt{h}*M_r*\eta_k]
      fhatr = @(x0,u,L) (Ehatr-h*Ahatr-h*Nhatr*u)\(x0+h*Bhatr*u+sqrt(h)*Mhatr*randn(size(Mhatr,2), L));
      disp(kk+" " + ii)
      rng(seed_test)
      
      % ROM simulation w/ testing I.C, input
      [Eopinf, Copinf, f1opinf, f2opinf] = compute_model(fhatr, Vr, ...
                                                            x0_test, u_test, L_test);
      
      % ROM vs FOM via relative error in E, C
      E_error_OpInf(ii,kk)  = norm(Eopinf - EFOM, "fro") / norm(EFOM, "fro");
      C_error_OpInf(ii,kk)  = page_norm(Copinf - CFOM) / page_norm(CFOM);
      % C_error_OpInf(ii,kk)  = norm(CROM - CFOM, 'fro') / norm(CFOM, 'fro');

      % weak error
      f1_error_OpInf(ii,kk) = abs(f1opinf - f1FOM) / abs(f1FOM);
      f2_error_OpInf(ii,kk) = abs(f2opinf - f2FOM) / abs(f2FOM);
    end
end



%% POD

E_error_POD = zeros(rmax, 1);
C_error_POD = zeros(rmax, 1);
f1_error_POD = zeros(rmax, 1);
f2_error_POD = zeros(rmax, 1);

for ii=1:rmax
    % disp(ii)

    Vr = Vrmax(:,1:ii);
    Er = Vr' * E * Vr;
    Ar = Vr' * A * Vr;
    Br = Vr' * B;
    Nr = Vr' * N * Vr;
    Mr = Vr' * M;
    fpod = @(x0,u,L) (Er-h*Ar-h*Nr*u)\(x0+h*Br*u + sqrt(h)*Mr*randn(size(Mr,2), L));

    disp("ROM dimension " + ii )
    rng(seed_test)

    [Epod,Cpod,f1pod,f2pod] = compute_model(fpod,Vr,x0_test,u_test,L_test);
    
    % Compute errors
    E_error_POD(ii)  = norm(Epod-EFOM, "fro") / norm(EFOM, "fro");
    C_error_POD(ii)  = page_norm(Cpod-CFOM) / page_norm(CFOM);
    f1_error_POD(ii) = abs(f1pod-f1FOM) / abs(f1FOM);
    f2_error_POD(ii) = abs(f2pod-f2FOM) / abs(f2FOM);
end


%% Compare Expectation (Eopinf)

% FOM expectation
[nx, Nt] = size(Eopinf);
x = linspace(0, 1, nx);          % spatial grid
t = linspace(0, (Nt-1)*h, Nt);   % time vector

figure()
hold on;
hLine1 = plot(x, EFOM(:,1), 'DisplayName', 'FOM');
hLine2 = plot(x, Eopinf(:,1), 'LineWidth', 1.5, 'DisplayName', 'OpInf ROM');
hLine3 = plot(x, Epod(:,1), 'LineWidth', 1.5, 'DisplayName', 'POD ROM');
xlabel('Space');
ylabel('Amplitude');
ylim([min(Eopinf(:)) max(Eopinf(:))]);
grid on;
legend()
title(sprintf('OpInf ROM Expectation 1D Heat equation t = %.3f', t(1)));

for k = 1:Nt
    set(hLine1, 'YData', EFOM(:,k), 'DisplayName', 'FOM');
    set(hLine2, 'YData', Eopinf(:,k), 'DisplayName', 'OpInf ROM');
    set(hLine3, 'YData', Epod(:,k), 'DisplayName', 'POD ROM');
    title(sprintf('OpInf ROM Expectation 1D Heat equation t = %.3f', t(k)));
    drawnow;
end

%% Compare covariance
Nt = size(Copinf, 3);
fps = 50;                 % desired frames per second
pause_time = 1/fps;       % time between frames

figure()
subplot(1,3,1)
im1 = imagesc(CFOM(:,:,1));
title('FOM Covariance');
axis equal tight
caxis([min(CFOM(:)) max(CFOM(:))]);
colormap(jet);
cb = colorbar();

subplot(1,3,2)
im2 = imagesc(Copinf(:,:,1));
title('OpInf ROM Covariance');
axis equal tight
caxis([min(Copinf(:)) max(Copinf(:))]);
colormap(jet);
cb = colorbar();

subplot(1,3,3)
im3 = imagesc(Cpod(:,:,1));
title('POD ROM Covariance');
axis equal tight
caxis([min(Cpod(:)) max(Cpod(:))]);
colormap(jet);
cb = colorbar();

for kk = 1:Nt

    set(im1, 'CData', CFOM(:,:,kk));
    set(im2, 'CData', Copinf(:,:,kk));
    set(im3, 'CData', Cpod(:,:,kk));

    subplot(1,3,1)
    title(sprintf('FOM'));
    subplot(1,3,2)
    title(sprintf('OpInf ROM'));
    subplot(1,3,3)
    title(sprintf('POD ROM'));

    sgtitle(sprintf('Covariance @ t = %.3f', t(kk)));
    drawnow;
    pause(pause_time);
end



%% save errors
error.E_error_OpInf = E_error_OpInf;
error.C_error_OpInf = C_error_OpInf;
error.f1_error_OpInf = f1_error_OpInf;
error.f2_error_OpInf = f2_error_OpInf;

error.E_error_POD = E_error_POD;
error.C_error_POD = C_error_POD;
error.f1_error_POD = f1_error_POD;
error.f2_error_POD = f2_error_POD;


%save( "./errors/error"+datestr(now, 'yyyy-mm-dd_HH-MM-SS')+".mat", '-struct', 'error', '-v7.3');

%% plot

ymin = 1e-6;

r = linspace(1, 10, 10);  % ROM dimensions

figure;

% Expectation
subplot(2,2,1)
plot(r, error.E_error_OpInf, '-o', 'LineWidth', 1.5, 'Color', [0 0.45 0.74]); hold on
plot(r, error.E_error_POD, '-s', 'LineWidth', 1.5, 'Color', [0.85 0.33 0.10]);
scatter(r, error.E_error_OpInf, 40, [0 0.45 0.74], 'filled');
scatter(r, error.E_error_POD, 40, [0.85 0.33 0.10], 'filled');
set(gca, 'YScale', 'log')
legend({'OpInf','POD'}, 'Location','best')
title("Expectation")
axis([1 rmax ymin 1])
grid on

% Covariance
subplot(2,2,2)
plot(r, error.C_error_OpInf, '-o', 'LineWidth', 1.5, 'Color', [0 0.45 0.74]); hold on
plot(r, error.C_error_POD, '-s', 'LineWidth', 1.5, 'Color', [0.85 0.33 0.10]);
scatter(r, error.C_error_OpInf, 40, [0 0.45 0.74], 'filled');
scatter(r, error.C_error_POD, 40, [0.85 0.33 0.10], 'filled');
set(gca, 'YScale', 'log')
legend({'OpInf','POD'}, 'Location','best')
title("Covariance")
axis([1 rmax ymin 1])
grid on

% f1
subplot(2,2,3)
plot(r, error.f1_error_OpInf, '-o', 'LineWidth', 1.5, 'Color', [0 0.45 0.74]); hold on
plot(r, error.f1_error_POD, '-s', 'LineWidth', 1.5, 'Color', [0.85 0.33 0.10]);
scatter(r, error.f1_error_OpInf, 40, [0 0.45 0.74], 'filled');
scatter(r, error.f1_error_POD, 40, [0.85 0.33 0.10], 'filled');
set(gca, 'YScale', 'log')
legend({'OpInf','POD'}, 'Location','best')
title("f1")
axis([1 rmax ymin 1])
grid on

% f2
subplot(2,2,4)
plot(r, error.f2_error_OpInf, '-o', 'LineWidth', 1.5, 'Color', [0 0.45 0.74]); hold on
plot(r, error.f2_error_POD, '-s', 'LineWidth', 1.5, 'Color', [0.85 0.33 0.10]);
scatter(r, error.f2_error_OpInf, 40, [0 0.45 0.74], 'filled');
scatter(r, error.f2_error_POD, 40, [0.85 0.33 0.10], 'filled');
set(gca, 'YScale', 'log')
legend({'OpInf','POD'}, 'Location','best')
title("f2")
axis([1 rmax ymin 1])
grid on


%
% norm(EFOM(:,end))
% norm(CFOM(:,:,end))



%%
norm(EFOM, "fro")
page_norm(CFOM)

norm(EROM, "fro")
page_norm(CROM)

%% Signal-to-noise (SNR) ratio (FOM states)
% SNR(t_i) = ||E^L_{f,i}||_2 / ||C^L_{f,i}||_F

EFOM_norm = vecnorm(EFOM, 2, 1);
CFOM_norm = squeeze(sum(CFOM.^2, [1 2]));
SNR = EFOM_norm' ./ CFOM_norm;

figure()
semilogy(EFOM_norm, 'linewidth', 2, 'linestyle', ':', ...
    'DisplayName', '$\|E_i^L\|_2$')
hold on
semilogy(CFOM_norm, 'linewidth', 2, 'linestyle', '--', ...
    'DisplayName', '$\|C_i^L\|_2$')
hold on
semilogy(SNR, 'linewidth', 2, 'DisplayName', 'SNR($t_i$)')
legend('Interpreter', 'latex', 'FontSize', 20)



%% Empirical moment snapshot matrix
%  F^L = [E^L, C^L] \in R^{n x (n+1)(s+1)}
E_L = EFOM;
C_L = reshape(CFOM, n, 100*1000);
F_L = horzcat(E_L, C_L);

[V2,S2,~] = svd(F_L, "econ");

%% Plot PODs computed in two ways
% X:   state snapshot
% F^L: moment snapshot
figure()

for ii=1:2:2*rmax
    row = (ii+1)/2;

    % left column: singular vectors
    subplot(rmax, 2, ii)
    plot(V(:,ii), 'DisplayName', "X snapshot")
    hold on
    plot(V2(:,ii), 'DisplayName', "$F^L$ snapshot")
    legend('Interpreter','latex')

    if row == 1
        title("Left singular vectors")
    end

    % right column: errors
    subplot(rmax, 2, ii+1)
    eigvec_sign_correction = sign(V(1,ii)/V2(1,ii));
    semilogy(abs(V(:,ii) - eigvec_sign_correction*V2(:,ii)))
    if row == 1
        title("Absolute errors")
    end
end


%%
function [y] = page_norm(A)
	y = sum(abs(A).^2, 'all');
	y = sqrt(y);
end
