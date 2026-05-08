clear; clc;

%% 0.25vpp

input = '0p25'; 

% Letters for the 10 files
letters = 'a':'j'; % 0p08, 0p15, 0p18, 0p25

% letters(ismember(letters, ['c','d','j'])) = []; % 0p04
% letters(ismember(letters, ['e'])) = []; % 0p07

% letters(ismember(letters, ['b', 'c', 'd', 'g'])) = []; % 0p08
% letters = ['b', 'c', 'd', 'g'];  % 0p08

% letters = 'a':'q';  % 0p10
% letters(ismember(letters, ['a', 'g'])) = []; 

% letters = 'a':'n'; % 0p20
% letters = 'a':'j';

% letters = 'a':'n'; % 0p30
% letters(ismember(letters, ['f','h', 'i', 'k', 'n'])) = [];

% letters = 'a':'p'; % 0p35
% letters(ismember(letters, ['m','n'])) = []; 


num_files = numel(letters);
base_path = sprintf('/disk/hyk049/DHM_new_1Dcenter/%s/', input);

first_file = sprintf('%sQ_1D_%svpp_b.h5', base_path, input);
t = h5read(first_file, '/t')';
x = h5read(first_file, '/x')';
nx = length(x);

% Preallocate cell array for all Q_1D signals
Q_1D = cell(1, num_files);

% Read all files
for k = 1:num_files
    filename = base_path + "Q_1D_" + input + "vpp_" + letters(k) + ".h5";
    Q_1D{k} = h5read(filename, '/Q_1D')';
end


% Split samples
split_size = 5760;
num_segs = floor(length(t)/split_size);
tt = t(1:split_size);

num_signals = numel(Q_1D);

Qstate_all = zeros(nx, num_segs*num_signals, split_size);
% Q_nom = cell(1, length(Q_1D));

for s = 1:num_signals
    % Truncate each experimental signal to an integer number of segments
    Q_experiment = Q_1D{s}(:, 1:split_size*num_segs);

    % Q_nom{s} = cell(1, num_segs);

    % Split each signal into segments and assign to Qstate_all
    for k = 1:num_segs
        idx_start = (k-1)*split_size + 1;
        idx_end   = k*split_size;

        % Compute column index range in Qstate_all
        col_idx = (s-1)*num_segs + k;

        Q_temp = Q_experiment(:, idx_start:idx_end);
        
        % Q_temp = Q_temp - Q_temp(:, 1);  % normalize each segment (start from zero)
        % Q_nom{s}{k} = Q_temp;

        % Assign the segment
        Qstate_all(:, col_idx, :) = Q_temp;
    end
end

% Plot expected dynamics
[T, X] = meshgrid(tt, x);

figure('Color','w');
surf(T, X, squeeze(mean(Qstate_all, 2)), 'EdgeColor', 'none');
xlabel('Time [s]'); ylabel('x [\mum]');
zlabel('Surface displacement [\mum]');
colorbar; colormap jet;


%%
exp = 10;

[T, X] = meshgrid(tt, x);
figure('Color','w', 'Position', [0 800 600 600]);
% surf(T, X, squeeze(mean(Qstate_all, 2)), 'EdgeColor', 'none');
surf(T, X, squeeze(mean(Qstate_all(:, num_segs*(exp-1)+1:num_segs*exp, :), 2)), 'EdgeColor', 'none');
% xlabel('Time [s]', 'FontSize', 30, 'Interpreter','latex');
% ylabel('x [$\mu$m]', 'FontSize', 30, 'Interpreter','latex');
% zlabel('Surface displacement [$\mu$m]', 'FontSize', 30, 'Interpreter','latex');
xlim([0 0.05]); xticks(0:0.01:0.05);
colormap jet; axis off;
ax = gca;
ax.FontSize = 22;
ax.TickLabelInterpreter = 'latex';
ax.XTickLabelRotation = 0;
% ax.ZTick = []; ax.XTick = []; ax.YTick = [];

%% Plot expected dynamics for each experiment

[T, X] = meshgrid(tt, x);

Z = cell(1, num_signals);

for s = 1:num_signals
    idx = (s-1)*num_segs + (1:num_segs);   % segments belonging to signal s
    Z{s} = squeeze(mean(Qstate_all(:, idx, :), 2));
end

% Compute common color limits
zmin = min(cellfun(@(z) min(z(:)), Z));
zmax = max(cellfun(@(z) max(z(:)), Z));

figure('Color','w', 'Position', [0 800 1900 600]);
tiledlayout(2, 5,'Padding','compact','TileSpacing','compact');

for k = 1:num_signals
    ax = nexttile;

    surf(ax, T, X, Z{k}, 'EdgeColor','none')

    xlabel(ax,'Time [s]')
    ylabel(ax,'x [\mum]')
    zlabel(ax,'Surface displacement [\mum]')

    xlim(ax,[0 0.06])
    xticks(ax,0:0.01:0.05)

    caxis(ax,[zmin zmax])
    view(ax,3)
    colormap(ax,jet)
    colorbar(ax)
end

%% FOM covariance
CFOM = page_cov(Qstate_all, true);


%% plot empirical covariance

Nt = size(CFOM, 3);
fps = 2000;
pause_time = 1/fps;

figure;
im = imagesc(CFOM(:,:,1));
colormap("jet");
colorbar;
caxis([min(CFOM(:)), max(CFOM(:))]);

for kk=1:Nt
    set(im, 'CData', CFOM(:,:,kk));
    title(sprintf('FOM Covariance t = %.2e', t(kk)));
    drawnow;
    pause(pause_time);
end


%% Filter signals

FS_DHM = 115200;
fc = 3e2;

% Design 4th-order Butterworth LPF
Wn = fc / (FS_DHM/2);
[b, a] = butter(4, Wn, 'low');

% Reshape to 2D for fast filtering
[nx, Nsamp, Nt] = size(Qstate_all);

Q_reshaped = reshape(Qstate_all, nx*Nsamp, Nt);

% Zero-phase filtering along time dimension
Q_filtered = filtfilt(b, a, Q_reshaped')';  % filtfilt works column-wise
% Q_filtered = filtfilt(b, a, Q_filtered')';  % filter twice for stability

% Reshape back to 3D
Qstate_filtered = reshape(Q_filtered, nx, Nsamp, Nt);

% Plot expected dynamics (filtered)
figure('color', 'w');
surf(T, X, squeeze(mean(Qstate_filtered, 2)), 'EdgeColor', 'none');
xlabel('Time [s]'); ylabel('x [\mum]');
zlabel('Surface displacement [\mum]');
colorbar; colormap jet;


%% Get POD basis

disp("Computing subspace Vr...")
[V,S,~] = svd(reshape(Qstate_all, nx, []), "econ");
disp("Subspace computation complete!")

% Stochastic OpInf
rmax = 20;

h = double(tt(2) - tt(1));
[~, L, s] = size(Qstate_all);

seed_test = 42;

% FOM expectation (average) over all samples
EFOM = squeeze(mean(Qstate_all, 2));
% EFOM = squeeze(mean(Qstate_filtered, 2));

% Testing noise sample numbers
L_test = L * [1];
% L_test = L * [1, 10, 100, 1000];



%%
t_scale = 0.1;

tt_test = tt*t_scale;
h_test = tt_test(2) - tt_test(1);

% Test with different number of noise samples

f0 = 7;

f_input = f0*1e6;
u_amp = 1;
u_train = u_amp*cos(2*pi*f_input*tt_test);


%% Test with different r

A_reg = 1;
N_reg = 1;
B_reg = 1;
% lambda = [A_reg, B_reg, N_reg];
lambda = [0, 0, 0];

% H_reg = logspace(0, 4, 20);
H_reg = 0;

isbilinear = true;

r_all = 5:20;

% Initialization
E_error = zeros(length(r_all), length(H_reg));
C_error = zeros(length(r_all), length(H_reg));

rng(seed_test)

for rr = 1:length(r_all)
    disp("r=" + r_all(rr))
    tic 

    Vr_temp = V(:, 1:r_all(rr));      % ROM basis
    xr0 = Vr_temp' * EFOM(:,1);       % ROM initial consdtion
    
    % Q_train_temp = pagemtimes(Vr_temp', Qstate_filtered);
    Q_train_temp = pagemtimes(Vr_temp', Qstate_all);       % R^{r x L x s}
    E_train = squeeze(mean(Q_train_temp, 2));              % R^{r x s}
    C_train = page_cov(Q_train_temp, true);                % R^{r x r x s}
    
    % Drift OpInf
    [Ehatr, Ahatr, Bhatr, Nhatr] = infer_drift_u(E_train, u_train, h_test, ...
                                                    isbilinear, lambda);

    for hh = 1:length(H_reg)

    % Diffusion OpInf
    [Mhatr, Khatr] = infer_diffusion_u(C_train, u_train, h_test, Ahatr, Nhatr, H_reg(hh));

    sigma = 1;

    % rng(seed_test)

    % noise = randn(size(Mhatr,2), L);
    noise_all = randn(size(Mhatr,2), L, s);
    % Noise_all{r_all(rr)} = noise_all;

    % % ROM function using Wiener process
    % fhatr = @(x0, u, L, k) ...
    %         (Ehatr - h_test*Ahatr - h_test*Nhatr*u) \ ...
    %         (x0 + h_test*Bhatr*u + sqrt(h_test)*Mhatr*sigma*noise_all(:,:,k)) ;

    % ROM function using Wiener process
    fhatr = @(x0, u, L, k, noise_k) ...
            (Ehatr - h_test*Ahatr - h_test*Nhatr*u) \ ...
            (x0 + h_test*Bhatr*u + sqrt(h_test)*Mhatr*sigma*noise_k) ;

   % Hurst = 0.2;
    % fhatr = @(x0, u, L) ...
    %         (Ehatr - h*Ahatr - h*Nhatr*u) \ ...
    %         (x0 + h*Bhatr*u + Mhatr*sigma*generate_fGn(size(Mhatr,2), L, h, Hurst)) ;
    
    % ROM simulation w/ testing I.C. across L samples
    [Eopinf, Copinf] = compute_model_u(fhatr, xr0, u_train, L, L, size(Mhatr,2));

    % Reconstruct FOM dimension
    E_recon = Vr_temp * Eopinf;  % R^{n x s}
    
    E_error(rr, hh) = norm(E_train - Eopinf, 'fro') / norm(E_train, 'fro');
    C_error(rr, hh) = page_norm(C_train - Copinf) / page_norm(C_train);

    end
    fprintf('Elapsed time: %.2f seconds.\n', toc); 

end

E_error_nonzero = E_error(E_error ~= 0);
[r_min, h_min] = find(E_error == min(E_error_nonzero));
disp("r optimal = " + r_all(r_min));
% disp("H_reg optimal = " + H_reg(h_min));

% plot r vs error
figure('Color','w');
plot(r_all, E_error, '-o', 'LineWidth', 2, ...
    'MarkerSize', 8, 'MarkerFaceColor', 'auto');
set(gca, 'YScale', 'log')
grid on; box on;
ylim([min(E_error), 1]);
xlabel('r (ROM dimension)', 'FontSize', 12)
ylabel('Error', 'FontSize', 12)
set(gca, 'FontSize', 12, 'LineWidth', 1.2)



%% Test with fixed u_amp, sigma + regularization

A_reg = logspace(-1, 3, 10);
B_reg = logspace(0, 8, 10);
% N_reg = logspace(0, 5, 10);
N_reg = logspace(-1, 5, 10);

B_reg = 0;

% H_reg_opt = H_reg(Hreg_min);
% H_reg_opt = 100;
H_reg = [1, 5, 10, 50, 100, 500, 1e3, 5e3, 1e4, 5e4];


% r_opt = 16;
r_opt = r_all(r_min);

Vr_temp = V(:, 1:r_opt);      % ROM basis
xr0 = Vr_temp' * EFOM(:,1);   % ROM initial consdtion

% Q_train_temp = pagemtimes(Vr_temp', Qstate_filtered);
Q_train_temp = pagemtimes(Vr_temp', Qstate_all);       % R^{r x L x s}
E_train = squeeze(mean(Q_train_temp, 2));              % R^{r x s}
C_train = page_cov(Q_train_temp, true);                % R^{r x r x s}

isbilinear = true;

if isbilinear
    outerLoop = 1:length(N_reg);
    ROM_form = 'non-autonomous bilinear ROM';
else
    outerLoop = 1:length(A_reg); 
    ROM_form = 'non-autonomous linear ROM';
end

% Initialization
Eopinfs = cell(1, length(N_reg));
EROM = cell(1, length(N_reg));
CROM = cell(1, length(N_reg));
E_error = zeros(length(N_reg), length(B_reg), length(H_reg));
C_error = zeros(length(N_reg), length(B_reg), length(H_reg));

H_norms = zeros(length(N_reg), length(B_reg), length(H_reg));

rng(seed_test)
for ii = outerLoop
    tic
    disp(ii + "th N_reg")

    Eopinfs{ii} = cell(1, length(B_reg));
    EROM{ii} = cell(1, length(B_reg));
    CROM{ii} = cell(1, length(B_reg));    
   
    for jj = 1:length(B_reg)

    Eopinfs{ii}{jj} = cell(1, length(H_reg));
    EROM{ii}{jj} = cell(1, length(H_reg));
    CROM{ii}{jj} = cell(1, length(H_reg));

    if isbilinear
        lambda = [0, B_reg(jj), N_reg(ii)];
    else
        lambda = [A_reg(ii), B_reg(jj), 0];
    end

    for kk=1:length(H_reg)
    H_reg_temp = H_reg(kk);

    % Drift OpInf
    [Ehatr, Ahatr, Bhatr, Nhatr] = infer_drift_u(E_train, u_train, h, ...
                                                    isbilinear, lambda); 
    % Diffusion OpInf
    [Mhatr, Khatr] = infer_diffusion_u(C_train, u_train, h, Ahatr, Nhatr, H_reg_temp);

    sigma = 1;

    noise_all = randn(size(Mhatr,2), L, s);

    % rng(seed_test)
    % ROM function using Wiener process
    fhatr = @(x0, u, L, k) ...
            (Ehatr - h*Ahatr - h*Nhatr*u) \ ...
            (x0 + h*Bhatr*u + sqrt(h)*Mhatr*sigma*noise_all(:,:,k)) ;

    
    % ROM simulation w/ testing I.C. across L samples
    % Eopinf: R^{r x s}, Copinf: R^{r x r x s} 
    [Eopinf, Copinf] = compute_model_u(fhatr, Vr_temp, xr0, u_train, L, L);

    Eopinfs{ii}{jj}{kk} = Eopinf;

    % Reconstruct FOM dimension
    E_recon = Vr_temp * Eopinf;  % R^{n x s}
     
    EROM{ii}{jj}{kk} = E_recon;
    CROM{ii}{jj}{kk} = Copinf;
    
    E_error(ii,jj,kk) = norm(E_train - Eopinf, 'fro') / norm(E_train, 'fro');
    C_error(ii,jj,kk) = page_norm(C_train - Copinf) / page_norm(C_train);

    H_norms(ii,jj,kk) = norm(Mhatr);
    
    end 

    end
    
    fprintf('Elapsed time: %.2f seconds.\n', toc); 
end 

%%
H_norm_max = squeeze(max(H_norms, [], 1));
E_error_min = squeeze(min(E_error, [], 1));

%%
% sigmas = [0.4, 1, 1.5, 3, 5, 10, 20, 30, 40, 50];
sigmas = [0.35, 0.5, 1, 2, 5, 10, 20, 30, 40, 55]; % 0p20

%%
figure('color', 'w', Position=[100, 900, 500, 600]);
subplot(2,1,1)
yyaxis left
plot(H_reg, H_norm_max, '-o');
set(gca, 'YScale', 'log');
ylabel('$\left\| \widetilde{M}_r \right\|_F$', 'Interpreter', 'latex');

yyaxis right
plot(H_reg, E_error_min, '--*');
% set(gca, 'YScale', 'log');
ylabel('$E_{\mathrm{error}}$', 'Interpreter', 'latex');

set(gca, 'XScale', 'log');
xlabel('$\lambda_{\widetilde{H}_r}$', 'Interpreter', 'latex');
grid on;

subplot(2,1,2)
plot(H_reg, sigmas, '--^');
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
xlabel('$\lambda_{\widetilde{H}_r}$', 'Interpreter', 'latex');
ylabel('noise adjustment $\sigma$', 'Interpreter', 'latex');
grid on;



%%
t_scale = 1e0;

tt_test = tt*t_scale;
h_test = tt_test(2) - tt_test(1);

% Test with different number of noise samples

f0 = 7;

f_input = f0*1e6;
u_amp = 1;
u_train = u_amp*cos(2*pi*f_input*tt_test);

%% Test with fixed u_amp, sigma + regularization

A_reg = logspace(-1, 3, 10);
B_reg = logspace(0, 8, 10);
N_reg = logspace(-1, 5, 10);

B_reg = 0;

H_reg_opt = 1e6;

sigma = 1;

r_opt = 16;
% r_opt = r_all(r_min);

Vr_temp = V(:, 1:r_opt);      % ROM basis
xr0 = Vr_temp' * EFOM(:,1);   % ROM initial consdtion

% Q_train_temp = pagemtimes(Vr_temp', Qstate_filtered);
Q_train_temp = pagemtimes(Vr_temp', Qstate_all);        % R^{r x L x s}
E_train = squeeze(mean(Q_train_temp, 2));               % R^{r x s}
C_train = page_cov(Q_train_temp, true);                 % R^{r x r x s}

isbilinear = true;

if isbilinear
    outerLoop = 1:length(N_reg);
    ROM_form = 'non-autonomous bilinear ROM';
else
    outerLoop = 1:length(A_reg); 
    ROM_form = 'non-autonomous linear ROM';
end

% Initialization
Eopinfs = cell(1, length(outerLoop));
EROM = cell(1, length(outerLoop));
CROM = cell(1, length(outerLoop));
E_error = zeros(length(outerLoop), length(B_reg));
C_error = zeros(length(outerLoop), length(B_reg));

rng(seed_test)


% Sample xr0 from the initial distribution N(E0,C0) instead of using a
% single point!!
% Fix initial condition sampling
C0 = C_train(:,:,1);
E0 = E_train(:,1);      % initial mean, r x 1

% Ensure PSD
C0 = (C0 + C0') / 2;
[V0, D0] = eig(C0);
d0_vec = max(real(diag(D0)), 0);

% Sample L initial conditions from N(E0, C0)
xr0_dist = real(V0) * diag(sqrt(d0_vec)) * randn(r_opt, L) + E0;  % r x L

% Large Monte Carlo sample from exact initial distribution
% L_large = 170;
% xr0_mc = real(V0) * diag(sqrt(d0_vec)) * randn(r_opt, L_large) + E0;

fprintf('xr0_dist cov norm:   %.2e\n', norm(cov(xr0_mc'), 'fro'));
fprintf('C_train(:,:,1) norm: %.2e\n', norm(C0, 'fro'));


% rng(seed_test)
for ii = outerLoop
    tic
    disp(ii + "th N_reg")

    E_opinfs = cell(1, length(B_reg));

    EROM{ii} = cell(1, length(B_reg));
    CROM{ii} = cell(1, length(B_reg));    
   
    for jj = 1:length(B_reg)

    if isbilinear
        lambda = [0, B_reg(jj), N_reg(ii)];
    else
        lambda = [A_reg(ii), B_reg(jj), 0];
    end

    % Drift OpInf
    [Ehatr, Ahatr, Bhatr, Nhatr] = infer_drift_u(E_train, u_train, h_test, ...
                                                    isbilinear, lambda); 
    % Diffusion OpInf
    [Mhatr, Khatr] = infer_diffusion_u(C_train, u_train, h_test, Ahatr, Nhatr, H_reg_opt);
    
    
fprintf('Mhatr norm: %.2e\n', norm(Mhatr));
fprintf('Hhat norm:  %.2e\n', norm(Mhatr*Mhatr')); 


    % noise_all = randn(size(Mhatr,2), L, s);

    % rng(seed_test)
    % ROM function using Wiener process

    fhatr = @(x0, u, L, k, noise_k) ...
        (Ehatr - h_test*Ahatr - h_test*Nhatr*u) \ ...
        (x0 + h_test*Bhatr*u + sqrt(h_test)*Mhatr*sigma*noise_k) ;

    % fhatr = @(x0, u, L, k) ...
    %         (Ehatr - h_test*Ahatr - h_test*Nhatr*u) \ ...
    %         (x0 + h_test*Bhatr*u + sqrt(h_test)*Mhatr*sigma*noise_all(:,:,k)) ;

    
    % ROM simulation w/ testing I.C. across L samples
    % Eopinf: R^{r x s}, Copinf: R^{r x r x s} 
    [Eopinf, Copinf, Xr] = compute_model_u(fhatr, xr0_dist, u_train, L_large, L_large, size(Mhatr,2));



% Time-resolved error
for k = [1, round(s/4), round(s/2), round(3*s/4), s]
    err = norm(C_train(:,:,k) - Copinf(:,:,k), 'fro') / norm(C_train(:,:,k), 'fro');
    fprintf('k=%4d | C_train: %.2e | Copinf: %.2e | rel err: %.2f\n', ...
        k, norm(C_train(:,:,k),'fro'), norm(Copinf(:,:,k),'fro'), err);
end


    Eopinfs{ii}{jj} = Eopinf;

    % Reconstruct FOM dimension
    E_recon = Vr_temp * Eopinf;  % R^{n x s}
     
    EROM{ii}{jj} = E_recon;
    CROM{ii}{jj} = Copinf;
    
    E_error(ii,jj) = norm(E_train - Eopinf, 'fro') / norm(E_train, 'fro');
    C_error(ii,jj) = page_norm(C_train - Copinf) / page_norm(C_train);
    
    end 
    
    fprintf('Elapsed time: %.2f seconds.\n', toc); 
end 

%
% Check reduced-state prediction
E_error_nonzero = E_error(E_error ~= 0);
[i_min, j_min] = find(E_error == min(E_error_nonzero));

% C_error_nonzero = C_error(C_error ~= 0);
% [i_min, j_min] = find(C_error == min(C_error_nonzero));

Eopinf_opt = Eopinfs{i_min}{j_min};

figure('Color','w', 'Position', [0 900 1400 600]); 
for r_check = 1:6
    subplot(2,3,r_check); hold on;

    plot(tt, E_train(r_check,:), 'LineWidth', 1.5, 'DisplayName','Experiment');
    plot(tt, Eopinf_opt(r_check,:), '--', 'LineWidth', 1.5, 'DisplayName','ROM');
    
    ylabel(['Coefficient r = ', num2str(r_check)], 'FontSize', 14);
    xlabel('Time [s]', 'FontSize', 14);
    % xlim([0,0.05]);
    ax = gca;
    ax.FontSize = 12;
    ax.LineWidth = 1.2;
    ax.Box = 'on';

    if r_check == 1
        legend('Location','best');
    end
    sgtitle(sprintf("Bilinear, H_{reg} = %d, \\sigma = %.1f", H_reg_opt, sigma), ...
        'FontSize', 16);
end   


%%
E_error_nonzero = E_error(E_error ~= 0);
[i_min, j_min] = find(E_error == min(E_error_nonzero));
EROM_opt = EROM{i_min}{j_min};
  
filename = sprintf('/disk/hyk049/WT_RomFit/%s/%svpp_r=%d_ABN_sigma%.1f_Hreg%d_f_%.1f.mat', ...
                   input, input, r_opt, sigma, H_reg_opt, f0)
%
save(filename, 'EFOM', 'EROM_opt');


%% Plot fitted result

[TT, X] = meshgrid(tt, x);

% Only consider nonzero entries
E_error_nonzero = E_error(E_error ~= 0);
[i_min, j_min] = find(E_error == min(E_error_nonzero));
EROM_opt = EROM{i_min}{j_min};

rel_error = abs(EFOM - EROM_opt) ./ max(abs(EFOM(:)));
% rel_error = abs(EFOM - EROM_opt) ./ abs(EFOM);
rel_error(:,1) = zeros(200,1);

figure('Color','w', 'Position', [0 800 1400 400]);
tile = tiledlayout(2, 3,'Padding','compact','TileSpacing','compact');

% Data containers
stateData = {EFOM, EROM_opt};
errorData = rel_error * 100;

% Row 1 (3D view)
ax = gobjects(1, 6);

ax(1) = nexttile; surf(TT,X,stateData{1},'EdgeColor','none');
title("Expected dynamics of"+newline+"experimental 1D capillary wave");

ax(2) = nexttile; surf(TT,X,stateData{2},'EdgeColor','none');
title(sprintf("Expected dynamics of"+newline+"stochastic " +  ROM_form));
% linkaxes(ax([1 2]), 'z')

ax(3) = nexttile; surf(TT,X,errorData,'EdgeColor','none');
title("Pointwise relative error (%)");

% Row 2 (Top view)
ax(4) = nexttile; surf(TT,X,stateData{1},'EdgeColor','none'); view(2);
ax(5) = nexttile; surf(TT,X,stateData{2},'EdgeColor','none'); view(2);
ax(6) = nexttile; surf(TT,X,errorData,'EdgeColor','none'); view(2);

% Formatting 
stateAxes = ax([1 2 4 5]);
errorAxes = ax([3 6]);

for a = stateAxes
    xlabel(a,"Time [s]"); ylabel(a,"x [\mum]");
    zlabel(a,"Surface elevation [\mum]");
    colormap(a,jet); colorbar(a);
end

for a = errorAxes
    xlabel(a,"Time [s]"); ylabel(a,"x [\mum]");
    zlabel(a,"Relative error (%)");
    colormap(a,hot); colorbar(a);
end

%  Shared color limits 
zmin = min([EFOM(:); EROM_opt(:)]);
zmax = max([EFOM(:); EROM_opt(:)]);

set(stateAxes,'CLim',[zmin zmax])
zlim(stateAxes,[zmin zmax])

sgtitle(sprintf("H_{reg} = %d, \\sigma = %.1f", H_reg_opt, sigma), ...
        'FontSize', 16);

% title(tile, sprintf("input: %s , $r=%d$", input, r_opt), 'Interpreter','latex');



%%
[T, X] = meshgrid(tt, x);

figure('Color','w', 'Position', [0 800 600 600]);
surf(T, X, EROM_opt, 'EdgeColor', 'none');
xlabel('Time [s]', 'FontSize', 20, 'Interpreter','latex');
ylabel('x [$\mu$m]', 'FontSize', 20, 'Interpreter','latex');
zlabel('Surface displacement [$\mu$m]', 'FontSize', 20, 'Interpreter','latex');
xlim([0 0.05]); xticks(0:0.01:0.05);
colormap jet;
ax = gca;
ax.FontSize = 16;
ax.TickLabelInterpreter = 'latex';
ax.XTickLabelRotation = 0;

%%
[T, X] = meshgrid(tt, x);

figure('Color','w', 'Position', [0 800 1500 500]);

% --- Subplot 1 (FOM)
ax1 = subplot(1,2,1);
surf(T, X, EFOM, 'EdgeColor', 'none');
hx = xlabel('Time [s]', 'FontSize', 30, 'Interpreter','latex');
hx.Rotation = 10;
hy = ylabel('Spatial dimension [$\mu$m]', 'FontSize', 20, 'Interpreter','latex');
hy.Rotation = -20;
zlabel('Surface displacement [$\mu$m]', 'FontSize', 20, 'Interpreter','latex');
title('Experiment (P=517.9mW)','Interpreter','latex');

% --- Subplot 2 (ROM)
ax2 = subplot(1,2,2);
surf(T, X, EROM_opt, 'EdgeColor', 'none');
hx = xlabel('Time [s]', 'FontSize', 30, 'Interpreter','latex');
hx.Rotation = 10;
hy = ylabel('Spatial dimension [$\mu$m]', 'FontSize', 20, 'Interpreter','latex');
hy.Rotation = -20;
% zlabel('Surface displacement [$\mu$m]', 'FontSize', 35, 'Interpreter','latex');
title('ROM prediction','Interpreter','latex');

% --- Shared z-axis limits
zmin = min([EFOM(:); EROM_opt(:)]);
zmax = max([EFOM(:); EROM_opt(:)]);
zlim(ax1, [zmin zmax]);
zlim(ax2, [zmin zmax]);

% --- Shared color scaling
caxis(ax1, [zmin zmax]);
caxis(ax2, [zmin zmax]);

% --- Link axes (optional but nice)
linkaxes([ax1, ax2], 'xyz');

% --- Colormap
colormap(jet);

% --- Single shared colorbar
cb = colorbar('Position',[0.92 0.15 0.02 0.7]);
cb.Label.String = 'Surface displacement [$\mu$m]';
cb.Label.Interpreter = 'latex';
cb.FontSize = 16;

% --- Formatting
set([ax1, ax2], 'FontSize', 25, 'TickLabelInterpreter','latex');
set([ax1, ax2], 'XLim', [0 0.05], 'XTick', 0:0.01:0.05);



%% Test with varing sigma's & H regs

N_reg = logspace(-1, 5, 10);
B_reg = 0;

H_reg = [1, 5, 10, 50, 100, 500, 1e3, 5e3, 1e4, 5e4];
sigmas = [0.35, 0.5, 1, 2, 5, 10, 20, 30, 40, 55]; % 0p20


r_opt = 11;
% r_opt = r_all(r_min);

Vr_temp = V(:, 1:r_opt);      % ROM basis
xr0 = Vr_temp' * EFOM(:,1);   % ROM initial consdtion

% Q_train_temp = pagemtimes(Vr_temp', Qstate_filtered);
Q_train_temp = pagemtimes(Vr_temp', Qstate_all);       % R^{r x L x s}
E_train = squeeze(mean(Q_train_temp, 2));              % R^{r x s}
C_train = page_cov(Q_train_temp, true);                % R^{r x r x s}

isbilinear = true;

if isbilinear
    outerLoop = 1:length(N_reg);
    ROM_form = 'non-autonomous bilinear ROM';
else
    outerLoop = 1:length(A_reg); 
    ROM_form = 'non-autonomous linear ROM';
end

% Initialization
Eopinfs = cell(1, length(N_reg));
EROM = cell(1, length(N_reg));
CROM = cell(1, length(N_reg));
E_error = zeros(length(N_reg), length(B_reg), length(H_reg));
C_error = zeros(length(N_reg), length(B_reg), length(H_reg));

H_norms = zeros(length(N_reg), length(B_reg), length(H_reg));

rng(seed_test)
for ii = outerLoop
    tic
    disp(ii + "th N_reg")

    Eopinfs{ii} = cell(1, length(B_reg));
    EROM{ii} = cell(1, length(B_reg));
    CROM{ii} = cell(1, length(B_reg));    
   
    for jj = 1:length(B_reg)

    Eopinfs{ii}{jj} = cell(1, length(H_reg));
    EROM{ii}{jj} = cell(1, length(H_reg));
    CROM{ii}{jj} = cell(1, length(H_reg));

    if isbilinear
        lambda = [0, B_reg(jj), N_reg(ii)];
    else
        lambda = [A_reg(ii), B_reg(jj), 0];
    end

    for kk=1:length(H_reg)
    
    H_reg_temp = H_reg(kk);
    sigma = sigmas(kk);

    % Drift OpInf
    [Ehatr, Ahatr, Bhatr, Nhatr] = infer_drift_u(E_train, u_train, h, ...
                                                    isbilinear, lambda); 
    % Diffusion OpInf
    [Mhatr, Khatr] = infer_diffusion_u(C_train, u_train, h, Ahatr, Nhatr, H_reg_temp);

    noise_all = randn(size(Mhatr,2), L, s);

    % rng(seed_test)
    % ROM function using Wiener process
    fhatr = @(x0, u, L, k) ...
            (Ehatr - h*Ahatr - h*Nhatr*u) \ ...
            (x0 + h*Bhatr*u + sqrt(h)*Mhatr*sigma*noise_all(:,:,k)) ;

    
    % ROM simulation w/ testing I.C. across L samples
    % Eopinf: R^{r x s}, Copinf: R^{r x r x s} 
    [Eopinf, Copinf] = compute_model_u(fhatr, Vr_temp, xr0, u_train, L, L);

    Eopinfs{ii}{jj}{kk} = Eopinf;

    % Reconstruct FOM dimension
    E_recon = Vr_temp * Eopinf;  % R^{n x s}
     
    EROM{ii}{jj}{kk} = E_recon;
    CROM{ii}{jj}{kk} = Copinf;
    
    E_error(ii,jj,kk) = norm(E_train - Eopinf, 'fro') / norm(E_train, 'fro');
    C_error(ii,jj,kk) = page_norm(C_train - Copinf) / page_norm(C_train);

    H_norms(ii,jj,kk) = norm(Mhatr);
    
    end 

    end
    
    fprintf('Elapsed time: %.2f seconds.\n', toc); 
end 



%% This shows the expectation error does not change much if we adjust noise level
% accordingly to the varying Hreg

H_norm_max = squeeze(max(H_norms, [], 1));
E_error_min = squeeze(min(E_error, [], 1));

figure('color', 'w', Position=[100, 900, 500, 600]);
subplot(2,1,1)
% yyaxis left
plot(H_reg, H_norm_max, '-o');
set(gca, 'YScale', 'log');
ylabel('$\left\| \widetilde{M}_r \right\|_F$', 'Interpreter', 'latex');

% yyaxis right
% plot(H_reg, E_error_min, '--*');
% % set(gca, 'YScale', 'log');
% ylabel('$E_{\mathrm{error}}$', 'Interpreter', 'latex');
% 
set(gca, 'XScale', 'log');
xlabel('$\lambda_{\widetilde{H}_r}$', 'Interpreter', 'latex');
grid on;

subplot(2,1,2)
yyaxis left
plot(H_reg, sigmas, '--^');
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
xlabel('$\lambda_{\widetilde{H}_r}$', 'Interpreter', 'latex');
ylabel('noise adjustment $\sigma$', 'Interpreter', 'latex');
grid on;

yyaxis right
plot(H_reg, E_error_min, '--*');
set(gca, 'YScale', 'log');
ylabel('$E_{\mathrm{error}}$', 'Interpreter', 'latex');



%% Plot covariance error progression
C_error_t = zeros(1,Nt);
for i=1:Nt
    C_error_t(i) = norm(C_train(:,:,i) - CROM{i_min}{j_min}(:,:,i), 'fro') ...
                        / norm(C_train(:,:,i), 'fro');
end

figure;
semilogy(tt, C_error_t);

%%
norm(C_train - CROM{i_min}{j_min}, 'fro') / norm(C_train, 'fro')



%% Test with corrected sigma, different r's

A_reg = logspace(-1, 3, 10);
B_reg = logspace(0, 8, 10);
N_reg = logspace(-1, 5, 10);

B_reg = 0;

% H_reg_opt = H_reg(Hreg_min);
H_reg_opt = 5760;


r_test = 5:20;

% Initialization
Eopinfs = cell(1, length(r_test));
EROM = cell(1, length(r_test));
CROM = cell(1, length(r_test));
E_error = cell(1, length(r_test));
C_error = cell(1, length(r_test));

% E_error = zeros(length(B_reg), length(N_reg));
% C_error = zeros(length(B_reg), length(N_reg));

rng(seed_test)
for i_r = 1:length(r_test)
    tic
    rr = r_test(i_r);
    disp("r = " + rr + " training..")

    Vr_temp = V(:, 1:rr);      % ROM basis
    xr0 = Vr_temp' * EFOM(:,1);   % ROM initial consdtion
    
    Q_train_temp = pagemtimes(Vr_temp', Qstate_all);       % R^{r x L x s}
    E_train = squeeze(mean(Q_train_temp, 2));              % R^{r x s}
    C_train = page_cov(Q_train_temp, true);                % R^{r x r x s}
    
    isbilinear = true;

    if isbilinear
        outerLoop = 1:length(N_reg);
        ROM_form = 'non-autonomous bilinear ROM';
    else
        outerLoop = 1:length(A_reg); 
        ROM_form = 'non-autonomous linear ROM';
    end

    % Initialization
    Eopinfs{i_r} = cell(1, length(N_reg));
    EROM{i_r} = cell(1, length(N_reg));
    CROM{i_r} = cell(1, length(N_reg));
    E_error{i_r} = zeros(length(B_reg), length(N_reg));
    C_error{i_r} = zeros(length(B_reg), length(N_reg));

    for ii = outerLoop

        Eopintfs{i_r}{ii} = cell(1, length(B_reg));
    
        EROM{i_r}{ii} = cell(1, length(B_reg));
        CROM{i_r}{ii} = cell(1, length(B_reg));    
       
        for jj = 1:length(B_reg)

        if isbilinear
            lambda = [0, B_reg(jj), N_reg(ii)];
        else
            lambda = [A_reg(ii), B_reg(jj), 0];
        end

        % Drift OpInf
        [Ehatr, Ahatr, Bhatr, Nhatr] = infer_drift_u(E_train, u_train, h, ...
                                                        isbilinear, lambda); 
        % Diffusion OpInf
        [Mhatr, Khatr] = infer_diffusion_u(C_train, u_train, h, Ahatr, Nhatr, H_reg_opt);
    
        sigma = 6;
    
        % ROM function using Wiener process
        fhatr = @(x0, u, L) ...
                (Ehatr - h*Ahatr - h*Nhatr*u) \ ...
                (x0 + h*Bhatr*u + sqrt(h)*Mhatr*sigma*randn(size(Mhatr,2), L)) ;
    
        % ROM simulation w/ testing I.C. across L samples
        % Eopinf: R^{r x s}, Copinf: R^{r x r x s} 
        [Eopinf, Copinf] = compute_model_u(fhatr, Vr_temp, xr0, u_train, L, L);

        Eopinfs{i_r}{ii}{jj} = Eopinf;

        % Reconstruct FOM dimension
        E_recon = Vr_temp * Eopinf;  % R^{n x s}
         
        EROM{i_r}{ii}{jj} = E_recon;
        CROM{i_r}{ii}{jj} = Copinf;
        
        E_error{i_r}(ii,jj) = norm(E_train - Eopinf, 'fro') / norm(E_train, 'fro');
        C_error{i_r}(ii,jj) = page_norm(C_train - Copinf) / page_norm(C_train);
        
        end 
    end

    fprintf('Elapsed time: %.2f seconds.\n', toc); 
end

%%
r_error_min = [];
EROM_opt_all = cell(1, length(r_test));

for i_r = 1:length(r_test)
    rr = r_test(i_r);
    E_error_nonzero = E_error{i_r}(E_error{i_r} ~= 0);
    min_error = min(E_error_nonzero);
    r_error_min = [r_error_min, min_error];
    [i_min, j_min] = find(E_error_nonzero == min_error);
    
    EROM_opt_all{i_r} = EROM{i_r}{i_min}{j_min};
   
end

%%
EROM_stack = cat(3, EROM_opt_all{:});
filename = sprintf('/disk/hyk049/WT_RomFit/%s/%svpp_sigma%d_r_all.mat', ...
                   input, input, sigma);
save(filename, 'EFOM', 'EROM_stack');

%%
figure;
plot(r_test, r_error_min, '-o', 'LineWidth', 2);


%% Plot fitted result

r_opt = 7;

[TT, X] = meshgrid(tt, x);

EROM_opt = EROM_opt_all{r_opt};

rel_error = abs(EFOM - EROM_opt) ./ max(abs(EFOM(:)));
% rel_error = abs(EFOM - EROM_opt) ./ abs(EFOM);
rel_error(:,1) = zeros(200,1);

figure('Color','w', 'Position', [0 800 1400 400]);
tile = tiledlayout(2, 3,'Padding','compact','TileSpacing','compact');

% Data containers
stateData = {EFOM, EROM_opt};
errorData = rel_error * 100;

% Row 1 (3D view)
ax = gobjects(1, 6);

ax(1) = nexttile; surf(TT,X,stateData{1},'EdgeColor','none');
title("Expected dynamics of"+newline+"experimental 1D capillary wave");

ax(2) = nexttile; surf(TT,X,stateData{2},'EdgeColor','none');
title(sprintf("Expected dynamics of"+newline+"stochastic " +  ROM_form));
% linkaxes(ax([1 2]), 'z')

ax(3) = nexttile; surf(TT,X,errorData,'EdgeColor','none');
title("Pointwise relative error (%)");

% Row 2 (Top view)
ax(4) = nexttile; surf(TT,X,stateData{1},'EdgeColor','none'); view(2);
ax(5) = nexttile; surf(TT,X,stateData{2},'EdgeColor','none'); view(2);
ax(6) = nexttile; surf(TT,X,errorData,'EdgeColor','none'); view(2);

% Formatting 
stateAxes = ax([1 2 4 5]);
errorAxes = ax([3 6]);

for a = stateAxes
    xlabel(a,"Time [s]"); ylabel(a,"x [\mum]");
    zlabel(a,"Surface elevation [\mum]");
    colormap(a,jet); colorbar(a);
end

for a = errorAxes
    xlabel(a,"Time [s]"); ylabel(a,"x [\mum]");
    zlabel(a,"Relative error (%)");
    colormap(a,hot); colorbar(a);
end

%  Shared color limits 
zmin = min([EFOM(:); EROM_opt(:)]);
zmax = max([EFOM(:); EROM_opt(:)]);

set(stateAxes,'CLim',[zmin zmax])
zlim(stateAxes,[zmin zmax])

% title(tile, sprintf("input: %s , $r=%d$", input, r_opt), 'Interpreter','latex');



%% Test with different input amp, noise amplitudes

A_reg = logspace(-1, 3, 10);
B_reg = logspace(0, 8, 10);
N_reg = logspace(0, 5, 10);

B_reg = 1;

r_opt = 15;
Vr_temp = V(:, 1:r_opt);

xr0 = Vr_temp' * EFOM(:,1);  % ROM initial condition

Q_train_temp = pagemtimes(Vr_temp', Qstate_all);  % R^{r x L x s}
E_train = squeeze(mean(Q_train_temp, 2));         % R^{r x s}
C_train = page_cov(Q_train_temp, true);           % R^{r x r x s}
    
isbilinear = true;

if isbilinear
    ROM_form = 'non-autonomous bilinear ROM';
    outerReg = N_reg;
else
    ROM_form = 'non-autonomous linear ROM';
    outerReg = A_reg;
end

f_input = 7e6;
u_amp = [1, 10, 100];

sigma = [1, 5, 10, 15, 20, 30];

E_error_all = cell(length(u_amp), 1);
C_error_all = cell(length(u_amp), 1);
E_opinf_all = cell(length(u_amp), 1);


for i_u = 1:length(u_amp)

u_train = u_amp(i_u)*cos(2*pi*f_input*tt);

E_error_all{i_u} = zeros(length(outerReg), length(B_reg), length(sigma));
C_error_all{i_u} = zeros(length(outerReg), length(B_reg), length(sigma));
E_opinf_all{i_u} = cell(length(outerReg), 1);

rng(seed_test)
tic
for ii = 1:length(outerReg)
    % rng(seed_test)
    disp(ii + "th outer reg")

    E_opinf_all{i_u}{ii} = cell(length(B_reg), 1);

    for jj = 1:length(B_reg)

        E_opinf_all{i_u}{ii}{jj} = zeros(r_opt, s, length(sigma));
    
        if isbilinear
            lambda = [1, B_reg(jj), N_reg(ii)];
        else
            lambda = [A_reg(ii), B_reg(jj), 1];
        end

        % Drift OpInf
        [Ehatr, Ahatr, Bhatr, Nhatr] = infer_drift_u(E_train, u_train, h, ...
                                                        isbilinear, lambda);
        % Diffusion OpInf
        [Mhatr, Khatr] = infer_diffusion_u(C_train, u_train, h, Ahatr, Nhatr);

        % Loop over noise amplitudes
        for i_sig = 1:length(sigma)
            
            fhatr = @(x0, u, L) ...
                    (Ehatr - h*Ahatr - h*Nhatr*u) \ ...
                    (x0 + h*Bhatr*u + sqrt(h)*Mhatr*sigma(i_sig)*randn(size(Mhatr,2), L)) ;
            
            [Eopinf, Copinf] = compute_model_u(fhatr, Vr_temp, xr0, u_train, L, L);
            
            E_error_all{i_u}(ii, jj, i_sig) = norm(E_train - Eopinf, 'fro') / norm(E_train, 'fro');
            % C_error_all{i_u}(ii,jj, i_sig) = page_norm(C_train - Copinf) / page_norm(C_train);
            
            E_opinf_all{i_u}{ii}{jj}(:,:, i_sig) = Eopinf;
        end 
    end
end
toc

end



%%
EROM_amp = cell(length(u_amp), 1);

for i_u = 1:length(u_amp)
    EROM_amp{i_u} = cell(length(sigma), 1);
    for i_sig = 1:length(sigma)
        E_error_sig = E_error_all{i_u}(:,:,i_sig);
        [i_min, j_min] = find(E_error_sig == min(E_error_sig(:)));
        E_opinf_min = E_opinf_all{i_u}{i_min}{j_min}(:,:,i_sig);
    
        Erecon = Vr_temp * E_opinf_min;
        EROM_amp{i_u}{i_sig} = Erecon;
    end
end


%% Save all ROMs
amps   = {'amp1', 'amp10', 'amp100'};
sigmas = {'sigma1', 'sigma5', 'sigma10', 'sigma15', 'sigma20', 'sigma30'};

for i = 1:length(amps)
    for j = 1:length(sigmas)
        EROM_all.(amps{i}).(sigmas{j}) = EROM_amp{i}{j};
    end
end


%%
EROM_opt = EROM_all.amp10.sigma10;

save(['/disk/hyk049/WT_RomFit/0p20/' ...
    '0p20vpp_r=18_ABN_amp10_sigma10_fc_1e3.mat'], 'EFOM', 'EROM_opt');


%% Plot all 
[TT, X] = meshgrid(tt, x);

EROM1  = EROM_all.amp10.sigma1;
EROM10 = EROM_all.amp10.sigma10;
EROM15 = EROM_all.amp10.sigma15;
EROM20 = EROM_all.amp10.sigma20;

rel1  = abs(EFOM - EROM1)  ./ abs(EFOM);
rel10 = abs(EFOM - EROM10) ./ abs(EFOM);
rel15 = abs(EFOM - EROM15) ./ abs(EFOM);
rel20 = abs(EFOM - EROM20) ./ abs(EFOM);

figure('Color','w', 'Position', [0 900 1700 1000]);
tile = tiledlayout(6, 4,'Padding','compact','TileSpacing','compact');

% Row1 FOM 3D
ax1 = nexttile(1);  surf(TT,X,EFOM,'EdgeColor','none');  title("Expected dynamics of" ...
                                + newline + "experimental 1D capillary wave");
% Row2 FOM top view
ax5 = nexttile(5);  surf(TT,X,EFOM,'EdgeColor','none');  view(2);

% Row3 ROM 3D
ax9 = nexttile(9);  surf(TT,X,EROM1,'EdgeColor','none'); 
title(ax9, sprintf('stochastic %s\n(\\sigma = 1)', ROM_form));

ax10 = nexttile(10); surf(TT,X,EROM10,'EdgeColor','none');
title(ax10, sprintf('stochastic %s\n(\\sigma = 10)', ROM_form));

ax11 = nexttile(11); surf(TT,X,EROM15,'EdgeColor','none');
title(ax11, sprintf('stochastic %s\n(\\sigma = 15)', ROM_form));

ax12 = nexttile(12); surf(TT,X,EROM20,'EdgeColor','none'); 
title(ax12, sprintf('stochastic %s\n(\\sigma = 20)', ROM_form));

% Row4 ROM top view
ax13 = nexttile(13);  surf(TT,X,EROM1,'EdgeColor','none'); view(2);
ax14 = nexttile(14);  surf(TT,X,EROM10,'EdgeColor','none'); view(2);
ax15 = nexttile(15);  surf(TT,X,EROM15,'EdgeColor','none'); view(2);
ax16 = nexttile(16);  surf(TT,X,EROM20,'EdgeColor','none'); view(2);

% Row5 error 3D
ax17 = nexttile(17);  surf(TT,X,rel1*100,'EdgeColor','none'); title("Pointwise relative" + ...
                                                            "error (%)");
ax18 = nexttile(18);  surf(TT,X,rel10*100,'EdgeColor','none'); title("Pointwise relative" + ...
                                                            "error (%)");
ax19 = nexttile(19);  surf(TT,X,rel15*100,'EdgeColor','none'); title("Pointwise relative" + ...
                                                            "error (%)");
ax20 = nexttile(20);  surf(TT,X,rel20*100,'EdgeColor','none'); title("Pointwise relative" + ...
                                                            "error (%)");
% Row6 error top view
ax21 = nexttile(21); surf(TT,X,rel1*100,'EdgeColor','none'); view(2);
ax22 = nexttile(22); surf(TT,X,rel10*100,'EdgeColor','none'); view(2);
ax23 = nexttile(23); surf(TT,X,rel15*100,'EdgeColor','none'); view(2);
ax24 = nexttile(24); surf(TT,X,rel20*100,'EdgeColor','none'); view(2);

% ===== Formatting (applied once) =====
allSurfAxes = [ax1 ax5 ax9 ax10 ax11 ax12 ax13 ax14 ax15 ax16];
allErrAxes  = [ax17 ax18 ax19 ax20 ax21 ax22 ax23 ax24];

for ax = allSurfAxes
    xlabel(ax,"Time [s]");
    ylabel(ax,"x [\mum]");
    colormap(ax,jet);
    colorbar(ax);
end

for ax = allErrAxes
    xlabel(ax,"Time [s]");
    ylabel(ax,"x [\mum]");
    colormap(ax,hot);
    colorbar(ax);
end

% Shared color limits for state plots
zmin = min([EFOM(:); EROM1(:); EROM10(:); EROM15(:); EROM20(:)]);
zmax = max([EFOM(:); EROM1(:); EROM10(:); EROM15(:); EROM20(:)]);

set(allSurfAxes,'CLim',[zmin zmax])
zlim(allSurfAxes,[zmin zmax])

% Shared error values
err_min = 0;
err_max = max([rel1(:); rel10(:); rel15(:); rel20(:)]) * 100;

set(allErrAxes,'CLim',[err_min err_max])
zlim(allErrAxes,[err_min err_max])

title(tile, sprintf("$r=%d$", r_opt), 'Interpreter','latex');



%% Test different input amplitude & norms of Ar, Br
f_input = 7e6;
u_amp = [1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3];
% u_train = u_amp*cos(2*pi*f_input*tt);

A_reg = logspace(-1, 3, 10);
B_reg = logspace(2, 6, 10);

A_norm_min = zeros(length(u_amp), 1);
B_norm_min = zeros(length(u_amp), 1);
EROM = cell(1, length(u_amp));
CROM = cell(1, length(u_amp));

E_error_min = zeros(length(u_amp), 1);

r_opt = 19;
Vr_temp = V(:, 1:r_opt);

Q_train_temp = pagemtimes(Vr_temp', Qstate_all);              % R^{r x L x s}
E_train = reshape(squeeze(mean(Q_train_temp, 2)), r_opt, s);  % R^{r x s}
C_train = page_cov(Q_train_temp, true);                       % R^{r x r x s}
    
isbilinear = false;

for i_u = 1:length(u_amp)

E_error = zeros(length(A_reg), length(B_reg));
C_error = zeros(length(A_reg), length(B_reg));
A_norm  = zeros(length(A_reg), length(B_reg));
B_norm  = zeros(length(A_reg), length(B_reg));

E_opinf_all = cell(1, length(A_reg));

tic
for ii = 1:length(A_reg)
    disp(ii + "th A_reg")  

    E_opinf_all{ii} = cell(1, length(B_reg));

    for jj = 1:length(B_reg)
    
    rng(seed_test)
    lambda = [A_reg(ii), B_reg(jj), 1];

    uu = u_amp(i_u)*cos(2*pi*f_input*tt);
    
    % Drift OpInf
    [Ehatr, Ahatr, Bhatr, Nhatr] = infer_drift_u(E_train, uu, h, isbilinear, lambda);
    
    % Diffusion OpInf
    [Mhatr, Khatr] = infer_diffusion_u(C_train, uu, h, Ahatr, Nhatr);

    A_norm(ii,jj) = norm(Ahatr);
    B_norm(ii,jj) = norm(Bhatr);
    
    sigma = 1;
    % ROM function using Wiener process
    fhatr = @(x0, u, L) ...
            (Ehatr - h*Ahatr - h*Nhatr*u) \ ...
            (x0 + h*Bhatr*u + sqrt(h)*Mhatr*sigma*randn(size(Mhatr,2), L)) ;
    
    % ROM initial consdtion
    xr0 = Vr_temp' * EFOM(:,1);
    
    % ROM simulation w/ testing I.C. across L samples
    % Eopinf: R^{r x s}, Copinf: R^{r x r x s} 
    [Eopinf, Copinf] = compute_model_u(fhatr, Vr_temp, xr0, u_train, L, L);
    
    E_opinf_all{ii}{jj} = Eopinf;
    
    E_error(ii,jj) = norm(E_train - Eopinf, 'fro') / norm(E_train, 'fro');
    C_error(ii,jj) = page_norm(C_train - Copinf) / page_norm(C_train);
    end
end
toc

[i_min, j_min] = find(E_error == min(E_error(:)));
EROM{i_u} = Vr_temp * E_opinf_all{i_min}{j_min};

E_error_min(i_u) = E_error(i_min, j_min); 
A_norm_min(i_u) = A_norm(9, 1);
B_norm_min(i_u) = B_norm(9, 1);

end

%%
figure('color', 'w'); hold on;

% Left y-axis
yyaxis left
h1 = scatter(u_amp, A_norm_min, 50, 'o');
hold on
h2 = scatter(u_amp, B_norm_min, 70, '*');
ylabel('Learned reduced operator norms')

% Right y-axis
yyaxis right
h3 = scatter(u_amp, E_error_min, '^');
ylabel('Expectation error')

% Log scales
set(gca, 'XScale', 'log')
yyaxis left
set(gca, 'YScale', 'log')
yyaxis right
set(gca, 'YScale', 'log')

xlabel('$\overline{u}$', 'Interpreter', 'latex', 'FontSize', 15)

lgd = legend([h1 h2 h3], ...
       {'$\|\widehat{A}_r\|_F$', '$\|\widehat{B}_r\|_F$', '$E_{error}$'}, ...
       'Interpreter', 'latex', ...
       'Location', 'best');

lgd.FontSize = 15;
grid on;


%%
save(['/disk/hyk049/WT_RomFit/0p25/' ...
    '0p25vpp_stoc_opinf_10datasets_B_10sigma.mat'], 'EFOM', 'EROM_opt');


%%
[TT, X] = meshgrid(tt, x);

% EROM_opt = EROM{9}{5};
EROM_opt = EROM{6};

figure('Color', 'w');
tile = tiledlayout(2, 3, 'Padding','compact','TileSpacing','compact');

% Subplot 1: FOM expectation
ax1 = nexttile;
s1 = surf(ax1, TT, X, EFOM, 'EdgeColor','none');
title(ax1,"Expected dynamics of" + newline + "experimental 1D capillary wave" ...
       + newline + "(L=170)");
xlabel(ax1,"Time [s]"); ylabel(ax1,"x [\mum]");
zlabel(ax1,"Surface elevation [\mum]");
% xlim(ax1,[0,0.06]); xticks(ax1,0:0.01:0.05);
colorbar(ax1); colormap(ax1, jet);

% Subplot 2: Stochastic OpInf ROM expectation
ax2 = nexttile;
s2 = surf(ax2, TT, X, EROM_opt, 'EdgeColor', 'none');
title(ax2,"Expected dynamics of" + newline + "stochastic linear nonautonomous ROM" ...
      + newline + "(L=170)");
xlabel(ax2,"Time [s]"); ylabel(ax2,"x [\mum]");
zlabel(ax2,"Surface elevation [\mum]");
% xlim(ax2,[0,0.06]); xticks(ax2,0:0.01:0.05);
colorbar(ax2); colormap(ax2, jet);


% Subplot 3: Pointwise relative error (%)
rel_error = abs(EFOM - EROM_opt)./abs(EFOM);

ax3 = nexttile;
s3 = surf(ax3, TT, X, rel_error*100, 'EdgeColor', 'none');
title(ax3,"Pointwise relative error (%)");
xlabel(ax3,"Time [s]"); ylabel(ax3,"x [\mum]");
zlabel(ax3,"Relative error (%)");
% xlim(ax3,[0,0.06]); xticks(ax3,0:0.01:0.05);
colorbar(ax3); colormap(ax3, hot);

ax4 = nexttile;
s4 = surf(ax4, TT, X, EFOM, 'EdgeColor','none');
view(ax4,2)
xlabel(ax4,"Time [s]"); ylabel(ax4,"x [\mum]");
zlabel(ax4,"Surface elevation [\mum]");
% xlim(ax4,[0,0.06]); xticks(ax4,0:0.01:0.05);
colorbar(ax4); colormap(ax4, jet);

ax5 = nexttile;
s5 = surf(ax5, TT, X, EROM_opt, 'EdgeColor', 'none');
view(ax5,2)
xlabel(ax5,"Time [s]"); ylabel(ax5,"x [\mum]");
zlabel(ax5,"Surface elevation [\mum]");
% xlim(ax5,[0,0.06]); xticks(ax5,0:0.01:0.05);
colorbar(ax5); colormap(ax5, jet);

ax6 = nexttile;
s6 = surf(ax6, TT, X, rel_error*100, 'EdgeColor', 'none');
view(ax6,2)
xlabel(ax6, "Time [s]"); ylabel(ax6,"x [\mum]");
zlabel(ax6,"Relative error (%)");
% xlim(ax6,[0,0.06]); xticks(ax6,0:0.01:0.05);
colorbar(ax6); colormap(ax6, hot);

title(tile, sprintf("$r=%d$", r_opt), 'interpreter', 'latex');

zmin = min([s1.ZData(:); s2.ZData(:); s4.ZData(:); s5.ZData(:)]);
zmax = max([s1.ZData(:); s2.ZData(:); s4.ZData(:); s5.ZData(:)]);
set([ax1 ax2 ax4 ax5], 'CLim', [zmin zmax])

zlim(ax1, [zmin zmax])
zlim(ax2, [zmin zmax])
zlim(ax4, [zmin zmax])
zlim(ax5, [zmin zmax])






%% 3 regularization

A_reg = logspace(-1, 3, 5);
B_reg = logspace(0, 8, 10);
N_reg = logspace(0, 5, 10);

r_opt = 14;
Vr_temp = V(:, 1:r_opt);

xr0 = Vr_temp' * EFOM(:,1);  % ROM initial condition

Q_train_temp = pagemtimes(Vr_temp', Qstate_all);  % R^{r x L x s}
E_train = squeeze(mean(Q_train_temp, 2));         % R^{r x s}
C_train = page_cov(Q_train_temp, true);           % R^{r x r x s}

isbilinear = true;

if isbilinear
    ROM_form = 'non-autonomous linear+bilinear ROM';
else
    ROM_form = 'non-autonomous linear ROM';
end

f_input = 7e6;
u_amp = [1, 10, 100];

sigma = [1, 5, 10, 15, 20, 30];

E_error_all = cell(length(u_amp), 1);
C_error_all = cell(length(u_amp), 1);
E_opinf_all = cell(length(u_amp), 1);

rng(seed_test)
for i_u = 1:length(u_amp)

u_train = u_amp(i_u)*cos(2*pi*f_input*tt);

E_error_all{i_u} = zeros(length(A_reg), length(B_reg), length(N_reg), length(sigma));
C_error_all{i_u} = zeros(length(A_reg), length(B_reg), length(N_reg), length(sigma));
E_opinf_all{i_u} = cell(length(A_reg), 1);

tic
for ii = 1:length(A_reg)
    disp(ii + "th A_reg-------------------")

    E_opinf_all{i_u}{ii} = cell(length(B_reg), 1);

    for jj = 1:length(B_reg)
        disp(jj + "th B_reg--------")
        E_opinf_all{i_u}{ii}{jj} = cell(length(N_reg), 1);

        for kk = 1:length(N_reg)
            disp(kk + "th N_reg")

            E_opinf_all{i_u}{ii}{jj}{kk} = zeros(r_opt, s, length(sigma));
    
            lambda = [A_reg(ii), B_reg(jj), N_reg(kk)];

            % Drift OpInf
            [Ehatr, Ahatr, Bhatr, Nhatr] = infer_drift_u(E_train, u_train, h, ...
                                                            isbilinear, lambda);
            % Diffusion OpInf
            [Mhatr, Khatr] = infer_diffusion_u(C_train, u_train, h, Ahatr, Nhatr);

            % Loop over noise amplitudes
            for i_sig = 1:length(sigma)
    
                fhatr = @(x0, u, L) ...
                        (Ehatr - h*Ahatr - h*Nhatr*u) \ ...
                        (x0 + h*Bhatr*u + sqrt(h)*Mhatr*sigma(i_sig)*randn(size(Mhatr,2), L)) ;
                
                [Eopinf, Copinf] = compute_model_u(fhatr, Vr_temp, xr0, u_train, L, L);
                
                E_error_all{i_u}(ii, jj, kk, i_sig) = norm(E_train - Eopinf, 'fro') / norm(E_train, 'fro');
                % C_error_all{i_u}(ii,jj, i_sig) = page_norm(C_train - Copinf) / page_norm(C_train);
                
                E_opinf_all{i_u}{ii}{jj}{kk}(:,:, i_sig) = Eopinf;
            end 
        end
    end
end
toc

end



%% Functions
function y = page_norm(A)
    y = sum(abs(A).^2, 'all');
    y = sqrt(y);
end

function fGn = generate_fGn(d, L, h, H)
    fGn = zeros(d, L);
    scale = h^H;
    for j = 1:d
        fBm = wfbm(H, L+1);   % unit time step
        fGn(j,:) = diff(fBm);
    end
    fGn = scale * fGn;        % correct scaling
end

function [stationary, drift] = decompose_LPF(x, fs, fc)
    % 3rd-order Butterworth low-pass filter
    [b, a] = butter(3, fc/(fs/2), 'low');
    
    % Zero-phase filtering (filtfilt)
    drift = filtfilt(b, a, x);
    
    % Stationary = original - drift
    stationary = x - drift;
end

