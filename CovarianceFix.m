clear; clc;

%% Load data

input = '0p35'; 

% Letters for the 10 files
% letters = 'a':'j'; % 0p08, 0p15, 0p18, 0p25
letters = 'a':'p';

% letters(ismember(letters, ['c','d','j'])) = []; % 0p04
% letters(ismember(letters, ['e'])) = []; % 0p07

% letters(ismember(letters, ['b', 'c', 'd', 'g'])) = []; % 0p08
% letters = ['b', 'c', 'd', 'g'];  % 0p08

% letters = 'a':'q';  % 0p10
% letters(ismember(letters, ['a', 'g'])) = []; 

% letters = 'a':'n'; % 0p20
% letters = 'a':'j';

% letters = 'a':'n'; % 0p30
% letters(ismember(letters, ['f','h', 'i', 'k' , 'n'])) = [];

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

    % Split each signal into segments and assign to Qstate_all
    for k = 1:num_segs
        idx_start = (k-1)*split_size + 1;
        idx_end   = k*split_size;

        % Compute column index range in Qstate_all
        col_idx = (s-1)*num_segs + k;

        Q_temp = Q_experiment(:, idx_start:idx_end);

        Q_temp = Q_temp;

        % Assign the segment
        Qstate_all(:, col_idx, :) = Q_temp;  
    end
end

%
% % Filter signals
% FS_DHM = 115200;
% fc = 1e4;
% 
% Qstate_all_low = LPF(Qstate_all, FS_DHM, fc);

% % Scale snapshots
% Qstate_all = Qstate_all/max(Qstate_all(:));

[T, X] = meshgrid(tt, x);
% Plot expected dynamics (filtered)
figure('color', 'w');
surf(T, X, squeeze(mean(Qstate_all, 2)), 'EdgeColor', 'none');
xlabel('Time [s]'); ylabel('x [\mum]');
zlabel('Surface displacement [\mum]');
colorbar; colormap jet;


% Get POD basis
disp("Computing subspace Vr...")
[V,S,~] = svd(reshape(Qstate_all, nx, []), "econ");
disp("Subspace computation complete!")

% Stochastic OpInf
rmax = 20;

h = double(tt(2) - tt(1));
[~, L, s] = size(Qstate_all);

% FOM expectation (average) over all samples
EFOM = squeeze(mean(Qstate_all, 2));

% Testing noise sample numbers
L_test = L * [1];



%%
t_scale = 1e-4;

tt_test = tt*t_scale;
h_test = tt_test(2) - tt_test(1);

% Test with different number of noise samples

f0 = 7;

f_input = f0*1e6;
u_amp = 1;
u_train = u_amp*cos(2*pi*f_input*tt_test);


%% Gaussianize
% STEP 0: Gaussianize Qstate_all before POD

[nx, L, s] = size(Qstate_all);

% Gaussianize each spatial DOF independently across L samples
% Qstate_all: nx x L x s
fprintf('Gaussianizing Qstate_all...\n');
Qstate_gauss = zeros(size(Qstate_all));
cdf_store_fom = cell(nx, 1);  % store CDF for each spatial DOF

for ii = 1:nx
    % Use samples at t=0 to define the transform
    row = squeeze(Qstate_all(ii, :, 1));  % 1 x L, samples at t=0
    [sorted, sort_idx] = sort(row);
    uniform = (1:L) / (L+1);
    
    % Store CDF for inversion later
    cdf_store_fom{ii}.sorted  = sorted;
    cdf_store_fom{ii}.uniform = uniform;
    
    % Apply same transform to ALL time steps
    for k = 1:s
        row_k = squeeze(Qstate_all(ii,:,k));  % 1 x L
        uniform_k = interp1(sorted, uniform, row_k, 'linear', 'extrap');
        uniform_k = max(min(uniform_k, 1-1e-6), 1e-6);  % clip to (0,1)
        Qstate_gauss(ii,:,k) = norminv(uniform_k);
    end
end

fprintf('Gaussianization complete!\n');
fprintf('Kurtosis before (spatial mean): %.2f\n', ...
    mean(abs(kurtosis(reshape(Qstate_all(:,:,1),  nx, L)'))));
fprintf('Kurtosis after  (spatial mean): %.2f\n', ...
    mean(abs(kurtosis(reshape(Qstate_gauss(:,:,1), nx, L)'))));

%% Get POD basis from Gaussianized data
disp("Computing subspace Vr from Gaussianized data...")
[V, S, ~] = svd(reshape(Qstate_gauss, nx, []), "econ");
disp("Subspace computation complete!")

% Stochastic OpInf
rmax = 20;
h    = double(tt(2) - tt(1));
[~, L, s] = size(Qstate_gauss);

% FOM expectation over Gaussianized samples
EFOM_gauss = squeeze(mean(Qstate_gauss, 2));

% Testing noise sample numbers
L_test = L * [1];



%% Test with different r

A_reg = 0;
B_reg = 0;
N_reg = 0;

H_reg = 1e6;

isbilinear = true;
isBu = false;

if isbilinear
    if isBu
        lambda = [A_reg, B_reg, N_reg];
    else
        lambda = [A_reg, N_reg];
    end
else
    if isBu
        lambda = [A_reg, B_reg];
    else
        lambda = [A_reg];
    end
end

r_all = 5:20;

% Initialization
E_error = zeros(length(r_all), length(H_reg));
C_error = zeros(length(r_all), length(H_reg));

% rng(seed_test)
% 
% % Sample xr0 from the initial distribution N(E0,C0), not using a single point!!
% % Fix initial condition sampling
% C0 = C_train(:,:,1);             % initial covaraince,  r x r
% E0 = E_train(: , 1);             % initial expectation, r x 1
% 
% C0 = (C0 + C0') / 2;
% R = chol(C0 + 1e-10*eye(r_opt));
% 
% % Sample L Monte Carlo IC samples from N(E0, C0)
% L_mc = L;
% xr0_dist = E0 + R' * randn(r_opt, L_mc);  % x = E0 + Rt*xi, where xi~N(0, I)
% 
% % This ensures Cov(xr0) = Rt*Cov(xi)*R = Rt*R = C0 

rng(seed_test)

for rr = 1:length(r_all)
    disp("r=" + r_all(rr))
    tic 

    Vr_temp = V(:, 1:r_all(rr));      % ROM basis
    % xr0 = Vr_temp' * EFOM(:,1);       % ROM initial consdtion
    
    Q_train_temp = pagemtimes(Vr_temp', Qstate_all);       % R^{r x L x s}
    E_train = squeeze(mean(Q_train_temp, 2));              % R^{r x s}
    C_train = page_cov(Q_train_temp, true);                % R^{r x r x s}
    
    % Drift OpInf
    [Ehatr, Ahatr, Bhatr, Nhatr] = infer_drift_u(E_train, u_train, h_test, ...
                                                    isbilinear, isBu, lambda);

    for hh = 1:length(H_reg)

    % Diffusion OpInf
    [Mhatr, Khatr] = infer_diffusion_u(C_train, u_train, h_test, Ahatr, Nhatr, H_reg(hh));

    sigma = 1;

    % noise = randn(size(Mhatr,2), L);
    noise_all = randn(size(Mhatr,2), L, s);
    % Noise_all{r_all(rr)} = noise_all;

    % ROM function using Wiener process
    % fhatr = @(x0, u, L, k) ...
    %         (Ehatr - h_test*Ahatr - h_test*Nhatr*u) \ ...
    %         (x0 + h_test*Bhatr*u + sqrt(h_test)*Mhatr*sigma*noise_all(:,:,k)) ;

    fhatr = @(x0, u, L, k, noise_k) ...
            (Ehatr - h_test*Ahatr - h_test*Nhatr*u) \ ...
            (x0 + h_test*Bhatr*u + sqrt(h_test)*Mhatr*sigma*noise_k) ;

    % ROM simulation w/ testing I.C. across L samples
    xr0_exp = squeeze(Q_train_temp(:,:,1));
    [Eopinf, Copinf] = compute_model_u(fhatr, xr0_exp, u_train, L, L, size(Mhatr,2));
    % [Eopinf, Copinf] = compute_model_u(fhatr, xr0, u_train, L, L);

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



%%
t_scale = 1;  % 8e-3 is the minimum (exactly at Nyquist)

tt_test = tt*t_scale;
% tt_test = linspace(0, tt(end), length(tt)/t_scale);
h_test = tt_test(2) - tt_test(1);

% input signal
f0      = 7;
f_input = f0 * 1e6;
u_amp   = 1;
u_train = u_amp*cos(2*pi*f_input*tt_test) + 0;

pts_per_cycle = 1 / (h_test * f_input);
fprintf('h_test        = %.4e s\n',  h_test)
fprintf('f_sample      = %.4e Hz\n', 1/h_test)
fprintf('pts/cycle     = %.4f\n',    pts_per_cycle)
fprintf('Nyquist limit = NEED at least 2 pts/cycle (f_sample >= %.1f MHz)\n', 2*f_input/1e6)



%% Test with fixed u_amp, sigma + regularization

A_reg = logspace(-1, 3, 10);
B_reg = logspace(0, 8, 10);
% N_reg = logspace(-1, 5, 10);

N_reg = logspace(0, 10, 10); % for better covariance estimate

A_reg_test = 0;
B_reg = 0;

H_reg_opt = 1e5;

sigma = 1;

r_opt = 11;
% r_opt = r_all(r_min);

isbilinear = true;
isBu = true;

Vr_temp = V(:, 1:r_opt);      % ROM basis
% xr0 = Vr_temp' * EFOM(:,1);   % ROM initial condition

Q_train_temp = pagemtimes(Vr_temp', Qstate_all);        % R^{r x L x s}
E_train = squeeze(mean(Q_train_temp, 2));               % R^{r x s}
C_train = page_cov(Q_train_temp, true);                 % R^{r x r x s}

f1_exp = mean(vecnorm(Q_train_temp(:,:,end).^2));
f2_exp = mean(Q_train_temp(:,:,end).^3.*exp(Q_train_temp(:,:,end)),'all');

if isbilinear
    outerLoop = 1:length(N_reg);
    ROM_form = 'non-autonomous bilinear ROM';
else
    outerLoop = 1:length(A_reg); 
    ROM_form = 'non-autonomous linear ROM';
end

% Set loop ranges based on flags
if isBu
    B_loop = 1:length(B_reg);
else
    B_loop = 1;  % only one iteration needed
    B_reg  = 0;  % force B_reg to zero
end

% Initialization
Eopinfs  = cell(1, length(outerLoop));
Copinfs  = cell(1, length(outerLoop));
EROM     = cell(1, length(outerLoop));
CROM     = cell(1, length(outerLoop));
ROMs_all = cell(1, length(outerLoop));
E_error  = zeros(length(outerLoop), length(B_loop));
C_error  = zeros(length(outerLoop), length(B_loop));
f1_error = zeros(length(outerLoop), length(B_loop));
f2_error = zeros(length(outerLoop), length(B_loop));

% Lambda assembly
if isbilinear && isBu
    get_lambda = @(ii,jj) [A_reg_test, B_reg(jj), N_reg(ii)];
elseif isbilinear && ~isBu
    get_lambda = @(ii,jj) [A_reg_test, N_reg(ii)];
elseif ~isbilinear && isBu
    get_lambda = @(ii,jj) [A_reg(ii), B_reg(jj), 0];
else  % ~isbilinear && ~isBu
    get_lambda = @(ii,jj) [A_reg(ii)];
end

% rng(seed_test)
for ii = outerLoop
    tic
    disp(ii + "th N_reg")

    EROM{ii}     = cell(1, length(B_loop));
    CROM{ii}     = cell(1, length(B_loop));
    ROMs_all{ii} = cell(1, length(B_loop));
   
    for jj = B_loop

    lambda = get_lambda(ii, jj);

    % Drift OpInf
    [Ehatr, Ahatr, Bhatr, Nhatr] = infer_drift_u(E_train, u_train, h_test, ...
                                                    isbilinear, isBu, lambda); 
    % Diffusion OpInf
    [Mhatr, Khatr] = infer_diffusion_u(C_train, u_train, h_test, Ahatr, Nhatr, H_reg_opt);
    
    % Save ROMs
    ROMs_all{ii}{jj} = struct('Ahat', Ahatr, 'Bhat', Bhatr, 'Nhat', Nhatr, ...
                            'Mhat', Mhatr, 'Khat', Khatr);

    rng(seed_test)  % same noise for all search cases
    % ROM function using Wiener process
    fhatr = @(x0, u, L, k, noise_k) ...
        (Ehatr - h_test*Ahatr - h_test*Nhatr*u) \ ...
        (x0 + h_test*Bhatr*u + sqrt(h_test)*Mhatr*sigma*noise_k) ;

    
    % ROM simulation w/ testing I.C. across L samples
    % Eopinf: R^{r x s}, Copinf: R^{r x r x s} 
    xr0_exp = squeeze(Q_train_temp(:,:,1));
    [Eopinf, Copinf, f1opinf, f2opinf, Xr] = compute_model_u(fhatr, xr0_exp, u_train, ...
                                                            L, L, size(Mhatr,2));

    % [Eopinf, Copinf, Xr] = compute_model_u(fhatr, xr0_dist, u_train, ...
    %                                             L_mc, L_mc, size(Mhatr,2));

    
% % Check Phi eigenvalues over time
% for k = [1, round(s/4), round(s/2), round(3*s/4), s]
%     Phi_k = (Ehatr - h_test*Ahatr - h_test*Nhatr*u_train(k)) \ eye(r_opt);
%     fprintf('k=%4d | u=%.3f | max|eig(Phi)|: %.8f\n', ...
%         k, u_train(k), max(abs(eig(Phi_k))));
% end

% Error over time
for k = [1, round(s/4), round(s/2), round(3*s/4), s]
    C_err = norm(C_train(:,:,k) - Copinf(:,:,k), 'fro') / norm(C_train(:,:,k), 'fro');
    fprintf('k=%4d | C_train: %.2e | Copinf: %.2e | rel err: %.2f\n', ...
        k, norm(C_train(:,:,k),'fro'), norm(Copinf(:,:,k),'fro'), C_err);
end


    Eopinfs{ii}{jj} = Eopinf;
    Copinfs{ii}{jj} = Copinf;

    % Reconstruct FOM dimension
    E_recon = Vr_temp * Eopinf;  % R^{n x s}
     
    EROM{ii}{jj} = E_recon;
    CROM{ii}{jj} = Copinf;
    
    E_error(ii,jj) = norm(E_train - Eopinf, 'fro') / norm(E_train, 'fro');
    C_error(ii,jj) = page_norm(C_train - Copinf) / page_norm(C_train);

    f1_error(ii,jj) = abs(f1_exp - f1opinf)/abs(f1_exp);
    f2_error(ii,jj) = abs(f2_exp - f2opinf)/abs(f2_exp);
    
    end 
    
    fprintf('Elapsed time: %.2f seconds.\n', toc); 
end


% Check reduced-state prediction
E_error_nonzero = E_error(E_error ~= 0);
[i_Emin, j_Emin] = find(E_error == min(E_error_nonzero));

Eopinf_opt = Eopinfs{i_Emin}{j_Emin};

figure('Color','w', 'Position', [0 900 1400 600]);
for r_check = 1:min(6, r_opt)
    subplot(2,3,r_check); hold on;

    plot(tt, E_train(r_check,:), 'LineWidth', 1.5, 'DisplayName','Experiment');
    plot(tt, Eopinf_opt(r_check,:), '--', 'LineWidth', 1.5, 'DisplayName','ROM');
    
    ylabel(['Coefficient r = ', num2str(r_check)], 'FontSize', 14);
    xlabel('Time [s]', 'FontSize', 14);
    xlim([tt(1),tt(end)]);
    ax = gca;
    ax.FontSize = 12; ax.LineWidth = 1.2;
    ax.Box = 'on';

    if r_check == 1
        legend('Location','best');
    end
    sgtitle(sprintf("Bilinear, H_{reg} = %d, \\sigma = %.1f", H_reg_opt, sigma), ...
        'FontSize', 16);
end   


% Plot covariance error over time
C_error_nonzero = C_error(C_error ~= 0);
[i_Cmin, j_Cmin] = find(C_error == min(C_error_nonzero));
Copinf_min = Copinfs{i_Cmin}{j_Cmin};
C_err_min = zeros(1,length(tt));

E_error_nonzero = E_error(E_error ~= 0);
[i_Emin, j_Emin] = find(E_error == min(E_error_nonzero));
Eopinf_min = Eopinfs{i_Emin}{j_Emin};
E_err_min = zeros(1,length(tt));
for is = 1:length(tt)
    C_err_min(is) = norm(C_train(:,:,is) - Copinf_min(:,:,is), 'fro') / norm(C_train(:,:,is), 'fro');
    E_err_min(is) = norm(E_train(:,is) - Eopinf_min(:,is), 'fro') / norm(E_train(:,is), 'fro');
end

%
figure('Color','w'); hold on;
plot(tt, C_err_min, 'LineWidth', 1.5, 'DisplayName', ...
            sprintf('Covariance error (avg=%.3f)', mean(C_err_min)));
plot(tt, E_err_min, '-.', 'LineWidth', 1.5, 'DisplayName', ...
            sprintf('Expectation error (avg=%.3f)', mean(E_err_min)));
yline(0.20, '--', 'DisplayName', '20%', 'Color', 'r', 'LineWidth', 1.5);
set(gca, 'YScale', 'log');
ylabel('Error', 'FontSize', 14);
xlabel('Time [s]', 'FontSize', 14);
legend('show');
box on;

%
% Get best result
Copinf_opt = Copinfs{i_Cmin}{j_Cmin};

% Frobenius norm over time
C_train_norm = zeros(1,s);
Copinf_norm  = zeros(1,s);
rel_err      = zeros(1,s);

for k = 1:s
    C_train_norm(k) = norm(C_train(:,:,k), 'fro');
    Copinf_norm(k)  = norm(Copinf_opt(:,:,k), 'fro');
    rel_err(k)      = norm(C_train(:,:,k) - Copinf_opt(:,:,k),'fro') / ...
                      norm(C_train(:,:,k), 'fro');
end

%
figure('Color','w', 'Position', [0 900 1200 800]);
% Panel 1: Covariance norm over time
subplot(3,1,1);
plot(tt, C_train_norm, 'b', 'LineWidth', 1.5, 'DisplayName', sprintf('Experiment (%d samples)', L));
hold on;
plot(tt, Copinf_norm,  'r--','LineWidth', 1.5, 'DisplayName', sprintf('ROM (%d samples)', L));
xlabel('Time [s]', 'FontSize', 14);
ylabel('$\|\mathbf{C}_r(t)\|_F$', 'FontSize', 14, 'Interpreter','latex');
title('Covariance norm over time', 'FontSize', 14);
legend('Location', 'best'); grid on;
ax = gca; ax.FontSize = 12; ax.Box = 'on';

% Panel 2: Relative error over time
subplot(3,1,2);
plot(tt, rel_err, 'k', 'LineWidth', 1.5);
yline(0,  'r--', 'LineWidth', 1.5, 'Label', 'Perfect match');
% yline(0.1, 'r--', 'LineWidth', 1.5, 'Label', '10% error');
yline(0.20, 'r--', 'LineWidth', 1.5, 'Label', '20% error');
% ylim([min(rel_err), max(rel_err)*1.2]);
xlabel('Time [s]', 'FontSize', 14);
ylabel(['$\|\mathbf{C}_r{_\textnormal{exp}}(t) - \mathbf{C}_r{_\textnormal{rom}}(t)\|_F /' ...
         '\|\mathbf{C}_r{_\textnormal{exp}}(t)\|_F$'], ...
           'FontSize', 18, 'Interpreter', 'latex');
yscale('log')
title('Covariance relative error over time', 'FontSize', 14);
grid on;
ax = gca; ax.FontSize = 12; ax.Box = 'on';

% Panel 3: Normalized covariance growth
C_train_growth = C_train_norm / C_train_norm(1);
Copinf_growth  = Copinf_norm  / Copinf_norm(1);

subplot(3,1,3);
plot(tt, C_train_growth, 'b',  'LineWidth', 1.5, 'DisplayName', sprintf('Experiment (%d samples)', L));
hold on;
plot(tt, Copinf_growth,  'r--','LineWidth', 1.5, 'DisplayName', sprintf('ROM (%d samples)', L));
yline(1.0, 'k--', 'LineWidth', 1, 'DisplayName', 'No growth');
xlabel('Time [s]', 'FontSize', 14);
ylabel('$\|\mathbf{C}_r(t) / \mathbf{C}_r(0)\|_F$', 'FontSize', 14, 'Interpreter','latex');
title('Normalized covariance growth', 'FontSize', 14);
legend('Location', 'best'); grid on;
ax = gca; ax.FontSize = 12; ax.Box = 'on';

sgtitle(sprintf('Covariance comparison: Experiment vs ROM (C_{error}=%.2f%%)', ...
    100*page_norm(C_train-Copinf_opt)/page_norm(C_train)), 'FontSize', 16);



%% Save results
filename = save_opt_result(E_error, EROM, EFOM, input, r_opt, sigma, H_reg_opt, isBu)



%% STEP 1: Pre-process — transform actual ICs to Gaussian
xr0_actual = squeeze(Q_train_temp(:,:,1));  % r x L

% Store transform parameters for each mode
mu_ic    = mean(xr0_actual, 2);   % r x 1
std_ic   = std(xr0_actual, 0, 2); % r x 1

% Quantile transform: x -> z (exactly Gaussian)
xr0_gauss = zeros(size(xr0_actual));
cdf_store = cell(r_opt, 1);  % store empirical CDF for inversion

for ii = 1:r_opt
    row = xr0_actual(ii,:);
    [sorted, sort_idx] = sort(row);
    uniform = (1:L) / (L+1);  % avoid 0 and 1
    
    % Store for inversion
    cdf_store{ii}.sorted   = sorted;
    cdf_store{ii}.uniform  = uniform;
    
    % Map to Gaussian
    [~, rank_idx] = sort(sort_idx);
    xr0_gauss(ii,:) = norminv(rank_idx / (L+1));
end

fprintf('Kurtosis before: %.2f\n', mean(abs(kurtosis(xr0_actual'))));
fprintf('Kurtosis after:  %.2f\n', mean(abs(kurtosis(xr0_gauss'))));

% STEP 2: Learn ROM in Gaussian space
% Transform ALL training data to Gaussian space
Q_gauss = zeros(size(Q_train_temp));  % r x L x s
for ii = 1:r_opt
    for ll = 1:L
        % Apply same marginal transform to all time steps
        row = squeeze(Q_train_temp(ii, ll, :))';
        % Use linear interpolation from stored CDF
        Q_gauss(ii, ll, :) = interp1(cdf_store{ii}.sorted, ...
                                      norminv(cdf_store{ii}.uniform), ...
                                      row, 'linear', 'extrap');
    end
end

% Compute statistics in Gaussian space
E_train_gauss = squeeze(mean(Q_gauss, 2));       % r x s
C_train_gauss = page_cov(Q_gauss, true);          % r x r x s

% Learn operators in Gaussian space
[Ehatr_g, Ahatr_g, Bhatr_g, Nhatr_g] = infer_drift_u(E_train_gauss, u_train, ...
                                                        h_test, isbilinear, lambda);
[Mhatr_g, Khatr_g] = infer_diffusion_u(C_train_gauss, u_train, h_test, ...
                                         Ahatr_g, Nhatr_g, H_reg_opt);

fprintf('Ahat_gauss max real eig: %.2e\n', max(real(eig(Ahatr_g))));

% STEP 3: Simulate ROM in Gaussian space
fhatr_g = @(x0, u, L, k, noise_k) ...
    (Ehatr_g - h_test*Ahatr_g - h_test*Nhatr_g*u) \ ...
    (x0 + h_test*Bhatr_g*u + sqrt(h_test)*Mhatr_g*sigma*noise_k);

[Eopinf_g, Copinf_g, f1opinf_g, f2opinf_g, Xr_g] = compute_model_u(fhatr_g, xr0_gauss, ...
                                                u_train, L, L, size(Mhatr_g,2));

% STEP 4: Post-process — transform back to original space
% Invert quantile transform for each sample at each time step
Xr_physical = zeros(size(Xr_g));  % r x L x s

for ii = 1:r_opt
    for k = 1:s
        z_samples = Xr_g(ii, :, k);  % Gaussian samples at time k

        % Gaussian → uniform → original space
        uniform_samples = normcdf(z_samples);
        Xr_physical(ii,:,k) = interp1(cdf_store{ii}.uniform, ...
                                        cdf_store{ii}.sorted, ...
                                        uniform_samples, 'linear', 'extrap');
    end
end

% Compute statistics in physical space
E_physical = squeeze(mean(Xr_physical, 2));     % r x s
C_physical = zeros(r_opt, r_opt, s);
for k = 1:s
    X = Xr_physical(:,:,k);
    mu = mean(X, 2);
    Xc = X - mu;
    C_physical(:,:,k) = (Xc * Xc') / (L-1);
end

% STEP 5: Compare
fprintf('=== Physical space comparison ===\n');
for k = [1, round(s/4), round(s/2), round(3*s/4), s]
    C_err = norm(C_train(:,:,k) - C_physical(:,:,k),'fro') / ...
            norm(C_train(:,:,k),'fro');
    E_err = norm(E_train(:,k)   - E_physical(:,k))         / ...
            norm(E_train(:,k));
    fprintf('k=%4d | E_err: %.3f | C_err: %.3f\n', k, E_err, C_err);
end
fprintf('Overall C_error: %.3f\n', ...
    page_norm(C_train - C_physical) / page_norm(C_train));


%% Compare original vs. Gaussinized ROM
% Initialization - add Gaussian versions
Eopinfs_g = cell(1, length(outerLoop));
Copinfs_g  = cell(1, length(outerLoop));
E_error_g  = zeros(length(outerLoop), length(B_reg));
C_error_g  = zeros(length(outerLoop), length(B_reg));
f1_error_g = zeros(length(outerLoop), length(B_reg));
f2_error_g = zeros(length(outerLoop), length(B_reg));


H_reg_test = 1e6;
r_opt = 10;

Vr_temp = V(:, 1:r_opt);
Q_train_temp = pagemtimes(Vr_temp', Qstate_all);        % R^{r x L x s}
E_train = squeeze(mean(Q_train_temp, 2));               % R^{r x s}
C_train = page_cov(Q_train_temp, true);                 % R^{r x r x s}

f1_exp = mean(vecnorm(Q_train_temp(:,:,end).^2));
f2_exp = mean(Q_train_temp(:,:,end).^3.*exp(Q_train_temp(:,:,end)),'all');

% STEP 1: Pre-process — compute once outside loop (transform doesn't depend on lambda)
xr0_actual = squeeze(Q_train_temp(:,:,1));  % r x L
xr0_gauss  = zeros(size(xr0_actual));
cdf_store  = cell(r_opt, 1);

for ii = 1:r_opt
    row = xr0_actual(ii,:);
    [sorted, sort_idx] = sort(row);
    uniform = (1:L) / (L+1);
    cdf_store{ii}.sorted  = sorted;
    cdf_store{ii}.uniform = uniform;
    [~, rank_idx] = sort(sort_idx);
    xr0_gauss(ii,:) = norminv(rank_idx / (L+1));
end

% Transform ALL training snapshots to Gaussian space — also once outside loop
Q_gauss = zeros(size(Q_train_temp));
for ii = 1:r_opt
    for ll = 1:L
        row = squeeze(Q_train_temp(ii,ll,:))';
        Q_gauss(ii,ll,:) = interp1(cdf_store{ii}.sorted, ...
                                    norminv(cdf_store{ii}.uniform), ...
                                    row, 'linear', 'extrap');
    end
end

% Compute Gaussian-space statistics — also once outside loop
E_train_gauss = squeeze(mean(Q_gauss, 2));
C_train_gauss = page_cov(Q_gauss, true);

fprintf('Kurtosis before: %.2f\n', mean(abs(kurtosis(xr0_actual'))));
fprintf('Kurtosis after:  %.2f\n', mean(abs(kurtosis(xr0_gauss'))));

% Grid search loop
rng(seed_test)
for ii = outerLoop
    tic
    disp(ii + "th N_reg")

    EROM{ii}    = cell(1, length(B_reg));
    CROM{ii}    = cell(1, length(B_reg));

    for jj = 1:length(B_reg)

        if isbilinear
            lambda = [A_reg_test, B_reg(jj), N_reg(ii)];
        else
            lambda = [A_reg(ii), B_reg(jj), 0];
        end

        % --- Original pipeline ---
        [Ehatr, Ahatr, Bhatr, Nhatr] = infer_drift_u(E_train, u_train, ...
                                                        h_test, isbilinear, lambda);
        [Mhatr, Khatr] = infer_diffusion_u(C_train, u_train, h_test, ...
                                             Ahatr, Nhatr, H_reg_test);

        fhatr = @(x0, u, L, k, noise_k) ...
            (Ehatr - h_test*Ahatr - h_test*Nhatr*u) \ ...
            (x0 + h_test*Bhatr*u + sqrt(h_test)*Mhatr*sigma*noise_k);

        xr0_exp = squeeze(Q_train_temp(:,:,1));
        [Eopinf, Copinf, f1opinf, f2opinf, Xr] = compute_model_u(fhatr, xr0_exp, ...
            u_train, L, L, size(Mhatr,2));

        Eopinfs{ii}{jj} = Eopinf;
        Copinfs{ii}{jj} = Copinf;
        EROM{ii}{jj}    = Vr_temp * Eopinf;
        CROM{ii}{jj}    = Copinf;

        E_error(ii,jj)  = norm(E_train - Eopinf,'fro') / norm(E_train,'fro');
        C_error(ii,jj)  = page_norm(C_train - Copinf) / page_norm(C_train);
        f1_error(ii,jj) = abs(f1_exp - f1opinf) / abs(f1_exp);
        f2_error(ii,jj) = abs(f2_exp - f2opinf) / abs(f2_exp);

        % --- Gaussianization pipeline ---
        % STEP 2: Learn operators in Gaussian space (same lambda)
        [Ehatr_g, Ahatr_g, Bhatr_g, Nhatr_g] = infer_drift_u(E_train_gauss, u_train, ...
                                                                h_test, isbilinear, lambda);
        [Mhatr_g, Khatr_g] = infer_diffusion_u(C_train_gauss, u_train, h_test, ...
                                                 Ahatr_g, Nhatr_g, H_reg_opt);

        fprintf('Ahat_gauss max real eig: %.2e\n', max(real(eig(Ahatr_g))));

        % STEP 3: Simulate in Gaussian space
        fhatr_g = @(x0, u, L, k, noise_k) ...
            (Ehatr_g - h_test*Ahatr_g - h_test*Nhatr_g*u) \ ...
            (x0 + h_test*Bhatr_g*u + sqrt(h_test)*Mhatr_g*sigma*noise_k);

        [Eopinf_g, Copinf_g, f1opinf_g, f2opinf_g, Xr_g] = compute_model_u(fhatr_g, ...
            xr0_gauss, u_train, L, L, size(Mhatr_g,2));

        % STEP 4: Post-process back to physical space
        Xr_physical = zeros(size(Xr_g));
        for mm = 1:r_opt
            for k = 1:s
                z_samples = Xr_g(mm,:,k);
                uniform_samples = normcdf(z_samples);
                Xr_physical(mm,:,k) = interp1(cdf_store{mm}.uniform, ...
                                               cdf_store{mm}.sorted, ...
                                               uniform_samples, 'linear', 'extrap');
            end
        end

        % Compute physical-space statistics
        E_physical = squeeze(mean(Xr_physical, 2));
        C_physical = zeros(r_opt, r_opt, s);
        for k = 1:s
            X  = Xr_physical(:,:,k);
            Xc = X - mean(X,2);
            C_physical(:,:,k) = (Xc*Xc') / (L-1);
        end

        % f1, f2 in physical space
        X_end_phys = Xr_physical(:,:,end);
        f1opinf_phys = mean(vecnorm(X_end_phys.^2));
        f2opinf_phys = mean(X_end_phys.^3 ./ exp(X_end_phys), 'all');

        % Store Gaussian results
        Eopinfs_g{ii}{jj} = E_physical;
        Copinfs_g{ii}{jj} = C_physical;

        E_error_g(ii,jj)  = norm(E_train - E_physical,'fro') / norm(E_train,'fro');
        C_error_g(ii,jj)  = page_norm(C_train - C_physical) / page_norm(C_train);
        f1_error_g(ii,jj) = abs(f1_exp - f1opinf_phys) / abs(f1_exp);
        f2_error_g(ii,jj) = abs(f2_exp - f2opinf_phys) / abs(f2_exp);

        % Print comparison
        fprintf('N_reg=%.2e | E_err: %.3f vs %.3f (gauss) | C_err: %.3f vs %.3f (gauss)\n', ...
            N_reg(ii), E_error(ii,jj), E_error_g(ii,jj), ...
            C_error(ii,jj), C_error_g(ii,jj));
    end

    fprintf('Elapsed time: %.2f seconds.\n', toc);
end

% Compare best results
[E_best,  ~]      = min(E_error(:));
[E_best_g, ~]     = min(E_error_g(:));
[C_best,  ~]      = min(C_error(:));
[C_best_g, ~]     = min(C_error_g(:));

fprintf('\n=== Best results comparison ===\n');
fprintf('Original:       E_error=%.3f | C_error=%.3f\n', E_best,   C_best);
fprintf('Gaussianized:   E_error=%.3f | C_error=%.3f\n', E_best_g, C_best_g);


%
% Find best results for both original and Gaussianized
C_error_nonzero   = C_error(C_error ~= 0);
C_error_g_nonzero = C_error_g(C_error_g ~= 0);

[~, idx_C]    = min(C_error_nonzero);
[~, idx_Cg]   = min(C_error_g_nonzero);

[i_Cmin,  j_Cmin]  = find(C_error   == min(C_error_nonzero));
[i_Cming, j_Cming] = find(C_error_g == min(C_error_g_nonzero));

Copinf_opt  = Copinfs{i_Cmin}{j_Cmin};      % best original
Copinf_opt_g = Copinfs_g{i_Cming}{j_Cming}; % best Gaussianized

fprintf('Best original:     N_reg=%.2e | C_error=%.3f\n', ...
    N_reg(i_Cmin),  C_error(i_Cmin,  j_Cmin));
fprintf('Best Gaussianized: N_reg=%.2e | C_error=%.3f\n', ...
    N_reg(i_Cming), C_error_g(i_Cming, j_Cming));

% Compute norms
C_train_norm  = zeros(1,s);
Copinf_norm   = zeros(1,s);
Copinf_norm_g = zeros(1,s);
rel_err       = zeros(1,s);
rel_err_g     = zeros(1,s);

for k = 1:s
    C_train_norm(k)  = norm(C_train(:,:,k),        'fro');
    Copinf_norm(k)   = norm(Copinf_opt(:,:,k),     'fro');
    Copinf_norm_g(k) = norm(Copinf_opt_g(:,:,k),   'fro');
    rel_err(k)       = norm(C_train(:,:,k) - Copinf_opt(:,:,k),   'fro') / C_train_norm(k);
    rel_err_g(k)     = norm(C_train(:,:,k) - Copinf_opt_g(:,:,k), 'fro') / C_train_norm(k);
end

% Normalized growth
C_train_growth  = C_train_norm  / C_train_norm(1);
Copinf_growth   = Copinf_norm   / Copinf_norm(1);
Copinf_growth_g = Copinf_norm_g / Copinf_norm_g(1);

% Plot
figure('Color','w', 'Position', [0 900 1200 800]);

% Panel 1: Covariance norm over time
subplot(3,1,1);
plot(tt, C_train_norm,  'b',   'LineWidth', 1.5, 'DisplayName', sprintf('Experiment (%d samples)', L));
hold on;
plot(tt, Copinf_norm,   'r--', 'LineWidth', 1.5, 'DisplayName', 'ROM');
plot(tt, Copinf_norm_g, 'm-.', 'LineWidth', 1.5, 'DisplayName', 'ROM Gauss');
xlabel('Time [s]', 'FontSize', 14);
ylabel('$\|\mathbf{C}_r(t)\|_F$', 'FontSize', 14, 'Interpreter', 'latex');
title('Covariance norm over time', 'FontSize', 14);
legend('Location', 'best'); grid on;
ax = gca; ax.FontSize = 12; ax.Box = 'on';

% Panel 2: Relative error over time
subplot(3,1,2);
plot(tt, rel_err,   'r--', 'LineWidth', 1.5, 'DisplayName', ...
    sprintf('ROM (avg=%.3f)',       mean(rel_err)));
hold on;
plot(tt, rel_err_g, 'm-.', 'LineWidth', 1.5, 'DisplayName', ...
    sprintf('ROM Gauss (avg=%.3f)', mean(rel_err_g)));
yline(0.20, 'k--', 'LineWidth', 1.5, 'DisplayName', '20% error');
xlabel('Time [s]', 'FontSize', 14);
ylabel(['$\|\mathbf{C}_{r,\textnormal{exp}} - \mathbf{C}_{r,\textnormal{rom}}\|_F /' ...
        '\|\mathbf{C}_{r,\textnormal{exp}}\|_F$'], 'FontSize', 14, 'Interpreter', 'latex');
yscale('log');
title('Covariance relative error over time', 'FontSize', 14);
legend('Location', 'best'); grid on;
ax = gca; ax.FontSize = 12; ax.Box = 'on';

% Panel 3: Normalized covariance growth
subplot(3,1,3);
plot(tt, C_train_growth,  'b',   'LineWidth', 1.5, 'DisplayName', sprintf('Experiment (%d samples)', L));
hold on;
plot(tt, Copinf_growth,   'r--', 'LineWidth', 1.5, 'DisplayName', 'ROM');
plot(tt, Copinf_growth_g, 'm-.', 'LineWidth', 1.5, 'DisplayName', 'ROM Gauss');
yline(1.0, 'k--', 'LineWidth', 1.0, 'DisplayName', 'No growth', 'HandleVisibility', 'off');
xlabel('Time [s]', 'FontSize', 14);
ylabel('$\|\mathbf{C}_r(t)\|_F / \|\mathbf{C}_r(0)\|_F$', 'FontSize', 14, 'Interpreter', 'latex');
title('Normalized covariance growth', 'FontSize', 14);
legend('Location', 'best'); grid on;
ax = gca; ax.FontSize = 12; ax.Box = 'on';

sgtitle(sprintf('Covariance: Original (%.2f%%) vs Gaussianized (%.2f%%)', ...
    100*page_norm(C_train-Copinf_opt)  / page_norm(C_train), ...
    100*page_norm(C_train-Copinf_opt_g)/ page_norm(C_train)), 'FontSize', 16);


%% Track max eigs of Phir
Phir = zeros(1,length(tt));
for kk = 1:length(tt)
    Phi_k = (Ehatr - h_test*Ahatr - h_test*Nhatr*u_train(kk)) \ eye(r_opt);
    Phir(kk) = max(abs(Phi_k(:)));
end

figure;
semilogy(tt, Phir, 'b--','LineWidth', 1.5); 


%% Quantify IC vs noise contribution to covariance

fhatr_test = @(x0, u, L, k, noise_k) ...
        (Ehatr - h_test*Ahatr - h_test*Nhatr*u) \ ...
        (x0 + h_test*Bhatr*u + sqrt(h_test)*Mhatr*10*noise_k) ;
% Run 1: actual ICs + diffusion (full model)
[~, Copinf_full, ~] = compute_model_u(fhatr_test, xr0_exp, u_train, ...
                                        L, L, size(Mhatr,2));

% Run 2: actual ICs + no diffusion (IC contribution only)
fhatr_nodiff = @(x0, u, L, k, noise_k) ...
    (Ehatr - h_test*Ahatr - h_test*Nhatr*u) \ ...
    (x0 + h_test*Bhatr*u);

[~, Copinf_nodiff, ~] = compute_model_u(fhatr_nodiff, xr0_exp, u_train, ...
                                          L, L, size(Mhatr,2));

% Run 3: zero IC (single point E0) + diffusion (noise contribution only)
xr0_det = E_train(:,1);  % single deterministic IC
[~, Copinf_noise, ~] = compute_model_u(fhatr_test, xr0_det, u_train, ...
                                         L, L, size(Mhatr,2));

% Global norms
fprintf('=== Global covariance contributions ===\n');
fprintf('C_train norm:        %.2e (100%%)\n', page_norm(C_train));
fprintf('Copinf_full norm:    %.2e (%.1f%%)\n', page_norm(Copinf_full), ...
    100*page_norm(Copinf_full)/page_norm(C_train));
fprintf('Copinf_nodiff norm:  %.2e (%.1f%%)\n', page_norm(Copinf_nodiff), ...
    100*page_norm(Copinf_nodiff)/page_norm(C_train));
fprintf('Copinf_noise norm:   %.2e (%.1f%%)\n', page_norm(Copinf_noise), ...
    100*page_norm(Copinf_noise)/page_norm(C_train));

% Time-resolved contributions
fprintf('\n=== Time-resolved contributions ===\n');
fprintf('%6s | %10s | %10s | %10s | %10s | %8s | %8s\n', ...
    'k', 'C_train', 'Full', 'IC only', 'Noise only', 'IC frac', 'Noise frac');
for k = [1, round(s/4), round(s/2), round(3*s/4), s]
    C_tr  = norm(C_train(:,:,k),       'fro');
    C_fu  = norm(Copinf_full(:,:,k),   'fro');
    C_ic  = norm(Copinf_nodiff(:,:,k), 'fro');
    C_no  = norm(Copinf_noise(:,:,k),  'fro');
    fprintf('%6d | %10.2e | %10.2e | %10.2e | %10.2e | %7.1f%% | %7.1f%%\n', ...
        k, C_tr, C_fu, C_ic, C_no, 100*C_ic/C_tr, 100*C_no/C_tr);
end

% Plot
figure('Color','w', 'Position', [0 900 1000 400]);
C_train_norm  = squeeze(vecnorm(reshape(C_train,       r_opt^2, []), 2, 1));
Copinf_norm   = squeeze(vecnorm(reshape(Copinf_full,   r_opt^2, []), 2, 1));
Cnodiff_norm  = squeeze(vecnorm(reshape(Copinf_nodiff, r_opt^2, []), 2, 1));
Cnoise_norm   = squeeze(vecnorm(reshape(Copinf_noise,  r_opt^2, []), 2, 1));

plot(tt, C_train_norm,  'k',  'LineWidth', 2,   'DisplayName', 'C\_train');
hold on;
plot(tt, Copinf_norm,   'b--','LineWidth', 1.5, 'DisplayName', 'Full ROM');
plot(tt, Cnodiff_norm,  'r--','LineWidth', 1.5, 'DisplayName', 'IC only');
plot(tt, Cnoise_norm,   'g--','LineWidth', 1.5, 'DisplayName', 'Noise only');
xlabel('Time [s]', 'FontSize', 14);
ylabel('||C||_F', 'FontSize', 14);
title('Covariance contribution: IC vs Noise', 'FontSize', 14);
legend('Location', 'best'); grid on;
ax = gca; ax.FontSize = 12; ax.Box = 'on';


%% A stability
Phieigs = eig(Ahatr);

figure('color', 'w', 'Position', [100, 100, 800, 320]);

scatter(real(Phieigs), imag(Phieigs), 80, 'filled');
hold on;

xline(0, '--', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1);
yline(0, '--', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1);
xline(0, '--', 'Color', 'red', 'LineWidth', 1);
xlabel('Real part');
ylabel('Imaginary part');
title('Eigenvalues of Ahatr');
grid on;
ax = gca;
ax.GridAlpha = 0.3;

hold off;


%% reduced- state sample trajectories, histogram
figure('Color', 'w');
i_coeff = 1;
data = squeeze(Q_train_temp(i_coeff, :, :));  % (N_samples, N_time)
% data = squeeze(Q_gauss(i_coeff, :, :));       % Gaussianized
terminal_points = data(:, 1);                 % initial points (first time step)
final_points    = data(:, end);               % final points (last time step)

% --- Settings ---
label_size = 15;
tick_size  = 15;
line_alpha = 0.5;
line_color = [0.20, 0.53, 0.63];
hist_color = [0.85, 0.35, 0.22];
hist_frac   = 0.22;
gap         = 0.015;
left_margin = 0.08;
traj_width  = 0.55;  % slightly narrower to fit right histogram

% --- Histogram (left) ---
hx = axes('Position', [left_margin, 0.13, hist_frac*traj_width, 0.78]);
histogram(hx, terminal_points, ...
    'Orientation', 'horizontal', ...
    'FaceColor',   hist_color, ...
    'FaceAlpha',   0.7, ...
    'EdgeColor',   'w', ...
    'LineWidth',   0.5);
hx.XDir     = 'reverse';
hx.FontSize = tick_size;
hx.XGrid    = 'on'; hx.YGrid = 'off';
hx.GridAlpha = 0.2; hx.GridLineWidth = 0.6;
xlabel(hx, 'Count', 'FontSize', label_size, 'Interpreter', 'latex');
title(hx, 'first time step', 'FontSize', label_size, 'FontWeight', 'normal', ...
             'Interpreter', 'latex');
box(hx, 'on');

% --- Trajectories (middle) ---
hx_pos = hx.Position;
ax = axes('Position', [hx_pos(1)+hx_pos(3)+gap, 0.13, traj_width, 0.78]);
hold(ax, 'on');
for i = 1:size(data, 1)
    plot(ax, tt, data(i,:), 'Color', [line_color, line_alpha], 'LineWidth', 0.6);
end
plot(ax, tt, mean(data, 1), 'Color', [0 0 0 1], 'linestyle', ':', 'LineWidth', 1.5);
xlabel(ax, 'Time [s]', 'FontSize', label_size, 'Interpreter', 'latex');
title(ax, sprintf('Coefficient %d (L = %d)', i_coeff, size(data,1)), ...
      'FontSize', 10, 'FontWeight','normal', 'Interpreter','latex');
ax.FontSize = tick_size;
yticks(ax, []);
grid(ax, 'on'); ax.GridAlpha = 0.2; ax.GridLineWidth = 0.6;
box(ax, 'on');

% --- Histogram (right) ---
ax_pos = ax.Position;
hx2 = axes('Position', [ax_pos(1)+ax_pos(3)+gap, 0.13, hist_frac*traj_width, 0.78]);
histogram(hx2, final_points, ...
    'Orientation', 'horizontal', ...
    'FaceColor',   hist_color, ...
    'FaceAlpha',   0.7, ...
    'EdgeColor',   'w', ...
    'LineWidth',   0.5);
hx2.FontSize = tick_size;
hx2.XGrid    = 'on'; hx2.YGrid = 'off';
hx2.GridAlpha = 0.2; hx2.GridLineWidth = 0.6;
hx2.YTickLabel = {};
xlabel(hx2, 'Count', 'FontSize', label_size, 'Interpreter', 'latex');
title(hx2, 'last time step', 'FontSize', label_size, 'FontWeight', 'normal', ...
        'Interpreter', 'latex');
box(hx2, 'on');

% --- Share y-axis across all three ---
ylim(hx,  ylim(ax));
ylim(hx2, ylim(ax));
hx.YTickLabel  = {};
ylabel(hx, '$\mathbf{X}_{r_1}$', 'FontSize', label_size, 'Interpreter', 'latex');
hx.YTickLabelMode  = 'auto';
hx.TickLabelInterpreter  = 'latex';
hx2.TickLabelInterpreter = 'latex';
ax.TickLabelInterpreter  = 'latex';


%% Save results
E_error_nonzero = E_error(E_error ~= 0);
[i_Emin, j_Emin] = find(E_error == min(E_error_nonzero));
EROM_opt = EROM{i_Emin}{j_Emin};

filename = sprintf('/disk/hyk049/WT_RomFit/%s/%svpp_r=%d_ABN_sigma%.1f_Hreg%d_f_%.1f.mat', ...
                   input, input, r_opt, sigma, H_reg_opt, f0)
%
save(filename, 'EFOM', 'EROM_opt');



%% Plot fitted result

[TT, X] = meshgrid(tt, x);

% Only consider nonzero entries
E_error_nonzero = E_error(E_error ~= 0);
[i_Emin, j_Emin] = find(E_error == min(E_error_nonzero));
EROM_opt = EROM{i_Emin}{j_Emin};

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


%%
t_scales = [8e-3, 1e-3, 1e-4, 1e-5, 1e-6];
N_reg    = logspace(0, 10, 10);
colors   = lines(length(t_scales));

A_reg_test = 0;
B_reg_test = 0;
H_reg_test = 1e5;

r_opt = 16;

Vr_temp = V(:, 1:r_opt);     
xr0 = Vr_temp' * EFOM(:,1); 

Q_train_temp = pagemtimes(Vr_temp', Qstate_all);        % R^{r x L x s}
E_train = squeeze(mean(Q_train_temp, 2));               % R^{r x s}
C_train = page_cov(Q_train_temp, true);                 % R^{r x r x s}

f1_exp = mean(vecnorm(Q_train_temp(:,:,end).^2));
f2_exp = mean(Q_train_temp(:,:,end).^3./exp(Q_train_temp(:,:,end)),'all');

figure('Color','w');
ax1 = subplot(1,2,1); hold(ax1,'on');
ax2 = subplot(1,2,2); hold(ax2,'on');

for ts_idx = 1:length(t_scales)
    t_scale = t_scales(ts_idx);
    tt_test = tt * t_scale;
    hh_test  = tt_test(2) - tt_test(1);
    u_train = u_amp * cos(2*pi*f_input*tt_test);
    pts_per_cycle = 1 / (hh_test * f_input);
    fprintf('\n--- t_scale=%.0e | h=%.2e s | %.2f pts/cycle ---\n', ...
             t_scale, hh_test, pts_per_cycle);

    E_error_ts = zeros(length(N_reg), 1);
    C_error_ts = zeros(length(N_reg), 1);
    Eopinfs_ts = cell(length(N_reg), 1);
    Copinfs_ts = cell(length(N_reg), 1);

    for ii = 1:length(N_reg)
        lambda = [A_reg_test, B_reg_test, N_reg(ii)];
        [Ehatr, Ahatr, Bhatr, Nhatr] = infer_drift_u(E_train, u_train, hh_test, ...
                                                       isbilinear, lambda);
        [Mhatr, ~] = infer_diffusion_u(C_train, u_train, hh_test, Ahatr, Nhatr, H_reg_test);
        fhatr = @(x0, u, L, k, noise_k) ...
            (Ehatr - hh_test*Ahatr - hh_test*Nhatr*u) \ ...
            (x0 + hh_test*Bhatr*u + sqrt(hh_test)*Mhatr*sigma*noise_k);

        xr0_exp = squeeze(Q_train_temp(:,:,1));
        [Eopinf, Copinf, ~] = compute_model_u(fhatr, xr0_exp, u_train, L, L, size(Mhatr,2));
        Eopinfs_ts{ii} = Eopinf;
        Copinfs_ts{ii} = Copinf;
        E_error_ts(ii) = norm(E_train - Eopinf, 'fro') / norm(E_train, 'fro');
        C_error_ts(ii) = page_norm(C_train - Copinf) / page_norm(C_train);
        fprintf('  N_reg=%.2e | E_err=%.4f | C_err=%.4f\n', ...
                 N_reg(ii), E_error_ts(ii), C_error_ts(ii));
    end

    E_nz = E_error_ts(E_error_ts ~= 0);
    C_nz = C_error_ts(C_error_ts ~= 0);
    i_Ebest = find(E_error_ts == min(E_nz), 1);
    i_Cbest = find(C_error_ts == min(C_nz), 1);
    Eopinf_best = Eopinfs_ts{i_Ebest};
    Copinf_best = Copinfs_ts{i_Cbest};
    fprintf('  >> best E: N_reg=%.2e (E_err=%.4f)\n', N_reg(i_Ebest), E_error_ts(i_Ebest));
    fprintf('  >> best C: N_reg=%.2e (C_err=%.4f)\n', N_reg(i_Cbest), C_error_ts(i_Cbest));

    C_err = zeros(1, length(tt));
    E_err = zeros(1, length(tt));
    for is = 1:length(tt)
        C_err(is) = norm(C_train(:,:,is) - Copinf_best(:,:,is), 'fro') ...
                  / norm(C_train(:,:,is), 'fro');
        E_err(is) = norm(E_train(:,is) - Eopinf_best(:,is), 'fro') ...
                  / norm(E_train(:,is), 'fro');
    end

    lbl = sprintf('t\\_scale=%.0e | %.1f pts/cycle | avg=%.3f', ...
                   t_scale, pts_per_cycle, mean(C_err));

    plot(ax1, 1:length(tt), E_err, '-.', 'LineWidth', 1.5, ...
                'Color', colors(ts_idx,:), 'HandleVisibility', 'off');
    plot(ax2, 1:length(tt), C_err, '-',  'LineWidth', 1.5, ...
                'Color', colors(ts_idx,:), 'DisplayName', lbl);
end

yline(ax2, 0.20, '--', 'Color', 'r', 'LineWidth', 2.5, 'DisplayName', '20% threshold');

% Format both axes
for ax = [ax1, ax2]
    set(ax, 'YScale', 'log');
    xlabel(ax, 'Time step',      'FontSize', 14);
    ylabel(ax, 'Relative error', 'FontSize', 14);
    ylim(ax, [1e-4, 1e0]);
    box(ax, 'on');
    grid(ax, 'on');
end

title(ax1, 'Mean error',       'FontSize', 14, 'FontWeight', 'normal');
title(ax2, 'Covariance error', 'FontSize', 14, 'FontWeight', 'normal');
legend(ax2, 'Location', 'best', 'FontSize', 10);



%% Test varying MC sample sizes when MC IC is used

L_larges = [170, 1700, 17000];
colors   = lines(length(L_larges));

figure('Color','w'); hold on;

for lL_idx = 1:length(L_larges)
    L_mc = L_larges(lL_idx);
    fprintf('\n--- L_large = %d ---\n', L_mc);

    % Sample L_large ICs from N(E0, C0)
    C0 = C_train(:,:,1);
    E0 = E_train(:,1);
    C0 = (C0 + C0') / 2;
    R  = chol(C0 + 1e-10*eye(r_opt));
    xr0_dist = E0 + R' * randn(r_opt, L_mc);

    % --- Regularization search ---
    E_error_ts = zeros(length(N_reg), 1);
    C_error_ts = zeros(length(N_reg), 1);
    Eopinfs_ts = cell(length(N_reg), 1);
    Copinfs_ts = cell(length(N_reg), 1);

    for ii = 1:length(N_reg)
        lambda = [0, B_reg, N_reg(ii)];

        [Ehatr, Ahatr, Bhatr, Nhatr] = infer_drift_u(E_train, u_train, h_test, ...
                                                       isbilinear, lambda);
        [Mhatr, ~] = infer_diffusion_u(C_train, u_train, h_test, Ahatr, Nhatr, H_reg_opt);

        fhatr = @(x0, u, L, k, noise_k) ...
            (Ehatr - h_test*Ahatr - h_test*Nhatr*u) \ ...
            (x0 + h_test*Bhatr*u + sqrt(h_test)*Mhatr*sigma*noise_k);

        [Eopinf, Copinf, ~] = compute_model_u(fhatr, xr0_dist, u_train, ...
                                               L_mc, L_mc, size(Mhatr,2));

        Eopinfs_ts{ii} = Eopinf;
        Copinfs_ts{ii} = Copinf;

        E_error_ts(ii) = norm(E_train - Eopinf, 'fro') / norm(E_train, 'fro');
        C_error_ts(ii) = page_norm(C_train - Copinf) / page_norm(C_train);

        fprintf('  N_reg=%.2e | E_err=%.4f | C_err=%.4f\n', ...
                 N_reg(ii), E_error_ts(ii), C_error_ts(ii));
    end

    % Pick best E and C separately
    E_nz = E_error_ts(E_error_ts ~= 0);
    C_nz = C_error_ts(C_error_ts ~= 0);
    i_Ebest = find(E_error_ts == min(E_nz), 1);
    i_Cbest = find(C_error_ts == min(C_nz), 1);

    Eopinf_best = Eopinfs_ts{i_Ebest};
    Copinf_best = Copinfs_ts{i_Cbest};

    fprintf('  >> best E: N_reg=%.2e (E_err=%.4f)\n', N_reg(i_Ebest), E_error_ts(i_Ebest));
    fprintf('  >> best C: N_reg=%.2e (C_err=%.4f)\n', N_reg(i_Cbest), C_error_ts(i_Cbest));

    % Error over time
    C_err = zeros(1, length(tt));
    E_err = zeros(1, length(tt));
    for is = 1:length(tt)
        C_err(is) = norm(C_train(:,:,is) - Copinf_best(:,:,is), 'fro') ...
                  / norm(C_train(:,:,is), 'fro');
        E_err(is) = norm(E_train(:,is) - Eopinf_best(:,is), 'fro') ...
                  / norm(E_train(:,is), 'fro');
    end

    lbl_C = sprintf('Cov | L=%d | avg=%.3f', L_mc, mean(C_err));
    lbl_E = sprintf('E | L=%d | avg=%.3f', L_mc, mean(E_err));
    subplot(1,2,1)
    plot(1:length(tt), E_err, '--', 'LineWidth', 1.5, 'Color', colors(lL_idx,:), ...
                    'DisplayName', lbl_E);
    subplot(1,2,2)
    plot(1:length(tt), C_err, '-',  'LineWidth', 1.5, 'Color', colors(lL_idx,:), ...
                    'DisplayName', lbl_C);
    yline(0.20, '--', 'Color', 'r', 'LineWidth', 2.5, 'DisplayName', '20% threshold');
end

yline(0.20, '--', 'Color', 'r', 'LineWidth', 2.5, 'DisplayName', '20% threshold');
ylim([1e-5, 1]);
set(gca, 'YScale', 'log');
ylabel('Relative error', 'FontSize', 14);
xlabel('Time step', 'FontSize', 14);
legend('Location', 'best', 'FontSize', 10);
% title('Covariance (—) and Expectation (-.) error vs. L\_large', 'FontSize', 14);
box on; grid on;

%%
Peigs = eig(Phi_k);

figure('color', 'w', 'Position', [100, 100, 800, 320]);

scatter(real(Peigs), imag(Peigs), 80, 'filled', 'MarkerFaceColor', [0, 0.5020, 0.5020]);
hold on;

xline(0, '--', 'Color', 'red', 'LineWidth', 1);
yline(0, '--', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1);

% Draw unit circle for reference
theta = linspace(0, 2*pi, 500);
plot(cos(theta), sin(theta), 'k--', 'LineWidth', 0.8);

xlabel('Real part');
ylabel('Imaginary part');
title('Eigenvalues of \Phi');
grid on;
ax = gca;
ax.GridAlpha = 0.3;

hold off;

%%
% Track cumulative Phi product eigenvalue
Phi_prod = eye(r_opt);
max_eig_cum = zeros(1,s);

u_input = u_train;
for k = 1:s
    Phi_k    = (Ehatr - h_test*Ahatr - h_test*Nhatr* u_input(k)) \ eye(r_opt);
    Phi_prod = Phi_k * Phi_prod;
    max_eig_cum(k) = max(abs(eig(Phi_prod)));
end

figure('Color','w');
subplot(2,1,1);
semilogy(tt, max_eig_cum, 'b', 'LineWidth', 1.5);
yline(1.0, 'r--', 'LineWidth', 1.5);
xlabel('Time [s]'); ylabel('max|eig(\Phi_{1:k})|');
title('Cumulative Phi eigenvalue growth');
grid on;

subplot(2,1,2);
plot(tt, u_input, 'b', 'LineWidth', 1.5);
xlabel('Time [s]'); ylabel('u(t)');
title('Input signal');
grid on;

%%
Copinf_opt = Copinfs{i_Emin}{j_Emin};

alpha = 1.96;

figure('Color','w', 'Position', [0 900 1400 600]);
for r_check = 1:6
    subplot(2,3,r_check); hold on;
    
    % Data mean and std
    E_data = E_train(r_check,:);
    std_data = squeeze(sqrt(abs(C_train(r_check, r_check, :))))';
    
    % ROM mean and std
    E_rom = Eopinf_opt(r_check,:);
    std_rom = squeeze(sqrt(abs(Copinf_opt(r_check, r_check, :))))';
    
    % if r_check == 1
    %     alpha = 0.01;
    % else
    %     alpha = 0.5;
    % end

    % Plot data confidence interval (1 sigma)
    fill([tt, fliplr(tt)], ...
         [E_data + alpha*std_data, fliplr(E_data - alpha*std_data)], ...
         'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none', ...
         'DisplayName', sprintf('±%g\\sigma exprmt', alpha));
    
    % Plot ROM confidence interval (1 sigma)
    fill([tt, fliplr(tt)], ...
         [E_rom + alpha*std_rom, fliplr(E_rom - alpha*std_rom)], ...
         'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none', ...
         'DisplayName', sprintf('±%g\\sigma ROM', alpha));
    
    % Plot means
    plot(tt, E_data, 'b',   'LineWidth', 1.5, 'HandleVisibility', 'off');
    plot(tt, E_rom,  'r--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    
    xlim([0, 0.05]);
    ylabel(['r = ', num2str(r_check)], 'FontSize', 14);
    xlabel('Time [s]', 'FontSize', 14);
    ax = gca;
    ax.FontSize = 12;
    ax.LineWidth = 1.2;
    ax.Box = 'on';
    legend('Location', 'best', 'FontSize', 10);
end


%% Functions
function y = page_norm(A)
    y = sum(abs(A).^2, 'all');
    y = sqrt(y);
end

function [stationary, drift] = decompose_LPF(x, fs, fc)
    % 3rd-order Butterworth low-pass filter
    [b, a] = butter(3, fc/(fs/2), 'low');
    
    % Zero-phase filtering (filtfilt)
    drift = filtfilt(b, a, x);
    
    % Stationary = original - drift
    stationary = x - drift;
end

function Q_filtered = LPF(x, Fs, fc)
    % Design 4th-order Butterworth LPF
    Wn = fc / (Fs/2);
    [b, a] = butter(4, Wn, 'low');
    
    % Reshape to 2D for fast filtering
    [nx, Nsamp, Nt] = size(x);
    
    Q_reshaped = reshape(x, nx*Nsamp, Nt);
    
    % Zero-phase filtering along time dimension
    Q_filtered = filtfilt(b, a, Q_reshaped')';  % filtfilt works column-wise
    
    % Reshape back to 3D
    Q_filtered = reshape(Q_filtered, nx, Nsamp, Nt);

end

% Post-processing function (for later use after simulation)
function x_physical = gauss_to_physical(z_samples, cdf_store_fom, nx, L)
% Transform Gaussianized samples back to physical space
% z_samples: nx x L
    x_physical = zeros(size(z_samples));
    for ii = 1:nx
        uniform_samples = normcdf(z_samples(ii,:));
        uniform_samples = max(min(uniform_samples, 1-1e-6), 1e-6);
        x_physical(ii,:) = interp1(cdf_store_fom{ii}.uniform, ...
                                    cdf_store_fom{ii}.sorted, ...
                                    uniform_samples, 'linear', 'extrap');
    end
end

function [filename] = save_opt_result(E_error, EROM, EFOM, ...
                                    input, r_opt, sigma, H_reg_opt, isBu)
%SAVE_OPT_RESULT Find optimal EROM entry and save results
%
% Inputs:
%   E_error    - matrix of errors
%   EROM       - cell array of models (same indexing as E_error)
%   input      - string (case name)
%   r_opt      - optimal rank
%   sigma      - noise level
%   H_reg_opt  - optimal regularization
%   f0         - parameter
%   EFOM       - full-order model
%
% Outputs:
%   EROM_opt   - selected optimal ROM
%   i_min,j_min- indices of minimum error
%   filename   - saved file path

    % Remove zeros
    E_error_nonzero = E_error(E_error ~= 0);

    if isempty(E_error_nonzero)
        error('E_error contains only zeros.');
    end

    % Find minimum (first occurrence)
    min_val = min(E_error_nonzero);
    idx = find(E_error == min_val, 1, 'first');
    [i_min, j_min] = ind2sub(size(E_error), idx);

    % Extract optimal ROM
    EROM_opt = EROM{i_min}{j_min};

    if isBu
        ROM_form = 'ABN';
    else
        ROM_form = 'AN';
    end

    % Create filename
    filename = sprintf('/disk/hyk049/WT_RomFit/%s/%svpp_r=%d_%s_sigma%.1f_Hreg%d.mat', ...
                       input, input, r_opt, ROM_form, sigma, H_reg_opt);

    % Save results
    save(filename, 'EFOM', 'EROM_opt');

end