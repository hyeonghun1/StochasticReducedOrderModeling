clc; clear;

%%
Q_c = h5read(['/disk/hyk049/DHM_new_experiment/' ...
    'Q_subspace_0p10vpp_c.h5'], '/Q_subspace');

t = h5read('/disk/hyk049/DHM_new_experiment/Q_subspace_0p18vpp.h5', '/t');
x = h5read('/disk/hyk049/DHM_new_experiment/Q_subspace_0p18vpp.h5', '/x');
y = h5read('/disk/hyk049/DHM_new_experiment/Q_subspace_0p18vpp.h5', '/y');

t = double(t); x = double(x); y = double(y);

Q_c = permute(Q_c, [3 2 1]);

Q_mean_c = squeeze(mean(Q_c, 2));
[nx, Nt] = size(Q_mean_c);

%%
Q_a = h5read(['/disk/hyk049/DHM_new_experiment/' ...
    'Q_subspace_0p10vpp_a.h5'], '/Q_subspace');
Q_a = permute(Q_a, [3 2 1]);
% Q_mean_a = squeeze(mean(Q_a, 2));

%%
Q_b = h5read(['/disk/hyk049/DHM_new_experiment/' ...
    'Q_subspace_0p10vpp_b.h5'], '/Q_subspace');
Q_b = permute(Q_b, [3 2 1]);

%%
Q_d = h5read(['/disk/hyk049/DHM_new_experiment/' ...
    'Q_subspace_0p10vpp_d.h5'], '/Q_subspace');
Q_d = permute(Q_d, [3 2 1]);

%%
Q_e = h5read(['/disk/hyk049/DHM_new_experiment/' ...
    'Q_subspace_0p10vpp_e.h5'], '/Q_subspace');
Q_e = permute(Q_e, [3 2 1]);


%% Split data into temporally independent datasets
split_size = 5760; % # of snapshots with
% one segment of temporal independence

num_segs = floor(size(Q_a, 3)/split_size);
tt = t(1:split_size);

Q_split_a = cell(1, num_segs);
Q_all_a = Q_a(:,:,1:split_size*num_segs);

Q_split_b = cell(1, num_segs);
Q_all_b = Q_b(:,:,1:split_size*num_segs);

Q_split_c = cell(1, num_segs);
Q_all_c = Q_c(:,:,1:split_size*num_segs);

Q_split_d = cell(1, num_segs);
Q_all_d = Q_d(:,:,1:split_size*num_segs);

Q_split_e = cell(1, num_segs);
Q_all_e = Q_e(:,:,1:split_size*num_segs);

for k=1:num_segs
    idx_start = (k-1)*split_size + 1;
    idx_end = k*split_size;

    Q_split_a{k} = Q_all_a(:, :, idx_start:idx_end);
    Q_split_b{k} = Q_all_b(:, :, idx_start:idx_end);
    Q_split_c{k} = Q_all_c(:, :, idx_start:idx_end);
    Q_split_d{k} = Q_all_d(:, :, idx_start:idx_end);
    Q_split_e{k} = Q_all_e(:, :, idx_start:idx_end);
end

clear Q_all_a Q_all_b Q_all_c Q_all_d Q_all_e 
clear Q_a Q_b Q_c Q_d Q_e


%% New state snapshots
Qstate_a = cat(2, Q_split_a{:});
Qstate_b = cat(2, Q_split_b{:});
Qstate_c = cat(2, Q_split_c{:});
Qstate_d = cat(2, Q_split_d{:});
Qstate_e = cat(2, Q_split_e{:});

clear Q_split_a Q_split_b Q_split_c Q_split_d Q_split_e


%%
Qstate_all = cat(2, Qstate_a, Qstate_b, Qstate_c, Qstate_d, Qstate_e);

clear Qstate_a Qstate_b Qstate_c Qstate_d Qstate_e

%%
% Qstate_all = permute(Qstate_all, [2 1 3]);

%%
[T, X] = meshgrid(tt, x);

figure;
surf(T, X, squeeze(mean(Qstate_all, 2)), 'EdgeColor', 'none');
xlabel('Time [s]'); ylabel('x [\mum]');
zlabel('Surface displacement [\mum]');
xlim([0 0.06]); xticks(0:0.01:0.05); 
colorbar; colormap jet;


%% FOM covariance
% CFOM = page_cov(Qstate_all, true);


%% plot empirical covariance

% Nt = size(CFOM, 3);
% fps = 2000;
% pause_time = 1/fps;
% 
% figure;
% im = imagesc(CFOM(:,:,1));
% colormap("jet");
% colorbar;
% caxis([min(CFOM(:)), max(CFOM(:))]);
% 
% for kk=1:Nt
%     set(im, 'CData', CFOM(:,:,kk));
%     title(sprintf('FOM Covariance t = %.2e', t(kk)));
%     drawnow;
%     pause(pause_time);
% end


%% Get POD basis

disp("Computing subspace Vr...")
[V,S,~] = svd(reshape(Qstate_all, nx, []), "econ");
disp("Subspace computation complete!")


%% plot singular values

normalized_sigma = diag(S)/max(diag(S));

figure;
semilogy(normalized_sigma);
hold on;
scatter(1:length(diag(S)), normalized_sigma, 'filled'); hold off;
xlabel('index');
ylabel('Singular values')
grid on;


%% Plot 1D POD modes

rmax = 10;
Vrmax = V(:, 1:rmax);

figure;
for ii = 1:rmax
    subplot(rmax, 1, ii)
    plot(V(:,ii)); hold on;
    ylabel(sprintf("Mode %d", ii));
end

%% Get projected data (reduced states)

rmax = 20;
Vrmax = V(:, 1:rmax);

% Project the FOM observations
Q_train = pagemtimes(Vrmax', Qstate_all);  % R^{r x L x s}


%% Projection error
proj_error = zeros(1, size(Qstate_all, 1));
diagS = diag(S);
sum_sigma = sum(diagS.^2);

for rr = 1:200
    proj_error(rr) = sqrt(sum(diagS(rr+1:end).^2) / sum_sigma) ;
end

figure('Color','w');
subplot(2,2,2)
semilogy(linspace(1,200,200), proj_error)
hold on
scatter(linspace(1,200,200), proj_error, '*')
yscale('log')
xlabel('Reduced dimension $r$', 'Interpreter','latex', 'FontSize', 15)
ylabel('Projection error $\rho_r$', 'Interpreter','latex', 'FontSize', 15)
grid on

subplot(2,2,1)
semilogy(normalized_sigma);
hold on;
scatter(1:length(diag(S)), normalized_sigma, 'filled'); hold off;
xlabel('Reduced dimension $r$', 'Interpreter','latex', 'FontSize', 15)
ylabel('Normalized singular values', 'Interpreter','latex', 'FontSize', 15)
grid on;

subplot(2,2,4)
semilogy(linspace(1,200,200), proj_error)
hold on
scatter(linspace(1,200,200), proj_error, '*')
yscale('log')
xlim([0,20])
xlabel('Reduced dimension $r$', 'Interpreter','latex', 'FontSize', 15)
ylabel('Projection error $\rho_r$', 'Interpreter','latex', 'FontSize', 15)
grid on

subplot(2,2,3)
semilogy(normalized_sigma);
hold on;
scatter(1:length(diag(S)), normalized_sigma, 'filled'); hold off;
xlim([0,20])
xlabel('Reduced dimension $r$', 'Interpreter','latex', 'FontSize', 15)
ylabel('Normalized singular values', 'Interpreter','latex', 'FontSize', 15)
grid on;



%% plot projections
QFOM_mean = squeeze(mean(Qstate_all, 2));

r_list = [1 3 5 10 20];

figure('Color','w');
tile = tiledlayout(2, 3, 'TileSpacing','compact', 'Padding','compact');

% Plot FOM
ax(1) = nexttile;
s(1) = surf(T, X, QFOM_mean, 'EdgeColor','none');
title('FOM')
xlabel('Time [s]'); ylabel('x [\mum]');
zlabel('Surface displacement [\mum]');
xlim([0 0.06]); xticks(0:0.01:0.05);
view(2)   % optional: top view (remove if you want 3D)

% Plot ROMs
for k = 1:numel(r_list)
    r = r_list(k);

    Vrtemp = Vrmax(:, 1:r);
    Q_proj = Vrtemp * (Vrtemp' * reshape(Qstate_all, nx, []));
    Q_proj = reshape(Q_proj, nx, size(Qstate_all,2), size(Qstate_all,3));

    QROM_mean = squeeze(mean(Q_proj, 2));

    ax(k+1) = nexttile;
    s(k+1) = surf(T, X, QROM_mean, 'EdgeColor','none');
    title(['Projection (r = ' num2str(r) ')'])
    xlabel('Time [s]'); ylabel('x [\mum]');
    xlim([0 0.06]); xticks(0:0.01:0.05);
    view(2)   % optional

    clear Q_proj QROM_mean
end

% Shared colormap and color limits
colormap(jet)
allC = vertcat(s.CData);
clim = [min(allC(:)), max(allC(:))];
set(ax, 'CLim', clim)

% One shared colorbar
cb = colorbar;
cb.Layout.Tile = 'east';
cb.Label.String = 'Surface displacement [\mum]';



%% OpInf

% max r value
rmax = 20;

s = size(Qstate_all, 3);
L = size(Qstate_all, 2);
h = tt(2) - tt(1);
isbilinear = false;

% FOM expectation
EFOM = squeeze(mean(Qstate_all, 2));

% relative weak errors
Q_T = Qstate_all(:,:,end);  % final state for all L noises (R^{r x L})
f1FOM = mean(vecnorm(Q_T.^2));        
f2FOM = mean(Q_T.^3./exp(Q_T), 'all');

% Initialization
EROM = cell(1,rmax); CROM = cell(1,rmax);
f1ROM = zeros(1,rmax); f2ROM = zeros(1,rmax);

E_error = cell(1,rmax); C_error = zeros(1,rmax);
f1_error = zeros(1,rmax); f2_error = zeros(1,rmax);


%% Test generating fractal Brownian motion
% rng default
H = 0.01;
l = 17000;

figure;
fBm_temp = wfbm(H, l, 'plot');


%% Train ROM

H = 0.50;  % Hurst exponent
fGn_all = cell(1, 20);

seed_test = 42;
doReg = true;
reg1 = logspace(-2, 5, 15);

for ii = 1:rmax

    disp("ROM dimension " + ii)
    rng(seed_test)

    Vr_temp = V(:, 1:ii);
    
    % Project the FOM observations
    Q_train_temp = pagemtimes(Vr_temp', Qstate_all);  % R^{r x L x s}
    
    % Estimate mean/covariance of the reduced observations
    % mean of reduced states across L realizations
    E_train = reshape(squeeze(mean(Q_train_temp, 2)), ii, s);
    
    % covariance of reduced states
    C_train = page_cov(Q_train_temp, true);
    
    for j=1:length(reg1)
    % Drift OpInf
    [Ehatr, Ahatr, Nhatr] = infer_drift(E_train, h, isbilinear, s, reg1(j));
    
    % Diffusion OpInf
    [Mhatr, Khatr] = infer_diffusion(C_train, h, Ahatr, Nhatr);

    % ROM initial consdtion
    xr0 = Vr_temp' * EFOM(:,1);

    d_temp = size(Mhatr,2);

    % % ROM function using Wiener process
    Wiener_noise = randn(d_temp, L);
    fhatr = @(x0, L) ...
        (Ehatr - h*Ahatr) \ (x0 + sqrt(h)*Mhatr*Wiener_noise) ;

    % Pre-generate fractal Gaussian noise
    % fGn_temp = zeros(d_temp, L);
    % if d_temp > 0
    %     for j = 1:d_temp
    %         fBm = wfbm(H, L+1);
    %         fGn_temp(j,:) = diff(fBm);
    %         fGn_all{ii} = fGn_temp;
    %     end
    % end
    % 
    % fhatr = @(x0, L) ...
    %     (Ehatr - h*Ahatr) \ (x0 + (h^H)*Mhatr*fGn_temp) ;

    % ROM simulation w/ testing I.C.
    [Eopinf, Copinf, f1opinf, f2opinf] = estimate(fhatr, Vr_temp, xr0, s, L);

    EROM_recon = Vr_temp*Eopinf;

    EROM{ii} = EROM_recon;   CROM{ii} = Copinf;
    f1ROM(ii) = f1opinf; f2ROM(ii) = f2opinf;

    E_err = norm(EFOM - EROM_recon, "fro") / norm(EFOM, "fro");
    E_error{ii} = [E_error{ii}, E_err];
    % C_error(ii) = page_norm(Copinf - CFOM) / page_norm(CFOM);
    % E_error(ii) = norm(E_train - Eopinf, "fro") / norm(E_train, "fro");
    C_error(ii) = page_norm(C_train - Copinf) / page_norm(C_train);
    f1_error(ii) = abs(f1opinf - f1FOM) / abs(f1FOM);
    f2_error(ii) = abs(f2opinf - f2FOM) / abs(f2FOM);
    end
end

error.E_error = E_error;
error.C_error = C_error;
error.f1_error = f1_error;
error.f2_error = f2_error;



%% plot errors
r = linspace(1, rmax, rmax);

figure;
% expectation
subplot(1,2,1)
plot(r, error.E_error, '-o', 'linewidth', 1.5); hold on;
scatter(r, error.E_error, 40, [0 0.45 0.74], 'filled');
set(gca, 'YScale', 'log')
title("Expectation error")
axis([1 rmax min(E_error)*0.9 max(E_error)*1.1])
grid on

% Covariance
subplot(1,2,2)
plot(r, error.C_error, '-o', 'LineWidth', 1.5); hold on;
scatter(r, error.C_error, 40, [0 0.45 0.74], 'filled');
set(gca, 'YScale', 'log')
title("Covariance error")
axis([1 rmax min(C_error)*0.9 max(C_error)*1.1])
grid on

% % f1
% subplot(2,2,3)
% plot(r, error.f1_error, '-o', 'LineWidth', 1.5); hold on;
% scatter(r, error.f1_error, 40, [0 0.45 0.74], 'filled');
% set(gca, 'YScale', 'log')
% title("f1")
% axis([1 rmax min(f1_error)*0.9 max(f1_error)*1.1])
% grid on
% 
% % f2
% subplot(2,2,4)
% plot(r, error.f2_error, '-o', 'LineWidth', 1.5); hold on;
% scatter(r, error.f2_error, 40, [0 0.45 0.74], 'filled');
% set(gca, 'YScale', 'log')
% title("f2")
% axis([1 rmax min(f2_error)*0.9 max(f2_error)*1.1])
% grid on


%%
[r_grid, reg_grid] = meshgrid((1:20), reg1);

E_error_array = cell2mat(E_error(:))';
figure('Color','w')
% scatter(repelem(1:20, 15), E_error_array(:));
surf(r_grid, reg_grid, E_error_array, 'EdgeColor','none')
% view(2)
set(gca, 'YScale', 'log'); set(gca, 'ZScale', 'log');
colormap(hot)
set(gca, 'ColorScale','log')
colorbar

xlabel('Reduced dimension r')
ylabel('$\lambda_1$', 'Interpreter','latex')
zlabel('Relative ROM error')


%%
[r_min, idx_min] = min(cellfun(@min, E_error))

%%
r_opt = 14;
[r_min, idx_min] = min(E_error{r_opt})


%% ROM simulation w/ optimal r & regularizer

rr = 14;
reg_idx = 12;

Vr_temp = V(:, 1:rr);
    
% Project the FOM observations
Q_train_temp = pagemtimes(Vr_temp', Qstate_all);  % R^{r x L x s}

% Estimate mean/covariance of the reduced observations
% mean of reduced states across L realizations
E_train = reshape(squeeze(mean(Q_train_temp, 2)), rr, s);

% covariance of reduced states
C_train = page_cov(Q_train_temp, true);

% Drift OpInf
[Ehatr, Ahatr, Nhatr] = infer_drift(E_train, h, isbilinear, s, reg1(reg_idx));

% Diffusion OpInf
[Mhatr, Khatr] = infer_diffusion(C_train, h, Ahatr, Nhatr);

% ROM initial consdtion
xr0 = Vr_temp' * EFOM(:,1);

d_temp = size(Mhatr,2);

% % ROM function using Wiener process
Wiener_noise = randn(d_temp, L);
fhatr = @(x0, L) ...
    (Ehatr - h*Ahatr) \ (x0 + sqrt(h)*Mhatr*Wiener_noise) ;

% Pre-generate fractal Gaussian noise
% fGn_temp = zeros(d_temp, L);
% if d_temp > 0
%     for j = 1:d_temp
%         fBm = wfbm(H, L+1);
%         fGn_temp(j,:) = diff(fBm);
%         fGn_all{ii} = fGn_temp;
%     end
% end
% 
% fhatr = @(x0, L) ...
%     (Ehatr - h*Ahatr) \ (x0 + (h^H)*Mhatr*fGn_temp) ;

% ROM simulation w/ testing I.C.
[Eopinf, Copinf, f1opinf, f2opinf] = estimate(fhatr, Vr_temp, xr0, s, L);

EROM_recon = Vr_temp*Eopinf;


%%
[TT, X] = meshgrid(tt, x);

% r_opt = 10;
% EROM_opt = V(:,1:r_opt)*EROM{r_opt}; 

% EROM_opt = EROM{r_opt};

ERPM_opt = EROM_recon;

figure('Color','w');
tile = tiledlayout(2,3,'Padding','compact','TileSpacing','compact');

% Subplot 1: FOM expectation
ax1 = nexttile;
s1 = surf(ax1, TT, X, EFOM, 'EdgeColor','none');
title(ax1,"FOM expectation");
xlabel(ax1,"Time [s]"); ylabel(ax1,"x [\mum]");
zlabel(ax1,"Surface elevation [\mum]");
xlim(ax1,[0,0.06]); xticks(ax1,0:0.01:0.05);
colorbar(ax1); colormap(ax1, jet);

% Subplot 2: Stochastic OpInf ROM expectation
ax2 = nexttile;
s2 = surf(ax2, TT, X, EROM_opt, 'EdgeColor', 'none');
title(ax2,"Stochastic OpInf ROM expectation");
xlabel(ax2,"Time [s]"); ylabel(ax2,"x [\mum]");
zlabel(ax2,"Surface elevation [\mum]");
xlim(ax2,[0,0.06]); xticks(ax2,0:0.01:0.05);
colorbar(ax2); colormap(ax2, jet);

% Match color limits for FOM and ROM
cmin = min([s1.ZData(:); s2.ZData(:)]);
cmax = max([s1.ZData(:); s2.ZData(:)]);
caxis(ax1,[cmin cmax]);
caxis(ax2,[cmin cmax]);

% Subplot 3: Pointwise relative error (%)
rel_error = abs(EFOM - EROM_opt)./EFOM;

ax3 = nexttile;
s3 = surf(ax3, TT, X, rel_error*100, 'EdgeColor', 'none');
title(ax3,"Pointwise relative error (%)");
xlabel(ax3,"Time [s]"); ylabel(ax3,"x [\mum]");
zlabel(ax3,"Relative error (%)");
xlim(ax3,[0,0.06]); xticks(ax3,0:0.01:0.05);
colorbar(ax3); colormap(ax3, hot);

ax4 = nexttile;
s4 = surf(ax4, TT, X, EFOM, 'EdgeColor','none');
view(ax4,2)
xlabel(ax4,"Time [s]"); ylabel(ax4,"x [\mum]");
zlabel(ax4,"Surface elevation [\mum]");
xlim(ax4,[0,0.06]); xticks(ax4,0:0.01:0.05);
colorbar(ax4); colormap(ax4, jet);

ax5 = nexttile;
s5 = surf(ax5, TT, X, EROM_opt, 'EdgeColor', 'none');
view(ax5,2)
xlabel(ax5,"Time [s]"); ylabel(ax5,"x [\mum]");
zlabel(ax5,"Surface elevation [\mum]");
xlim(ax5,[0,0.06]); xticks(ax5,0:0.01:0.05);
colorbar(ax5); colormap(ax5, jet);

ax6 = nexttile;
s6 = surf(ax6, TT, X, rel_error*100, 'EdgeColor', 'none');
view(ax6,2)
xlabel(ax6, "Time [s]"); ylabel(ax6,"x [\mum]");
zlabel(ax6,"Relative error (%)");
xlim(ax6,[0,0.06]); xticks(ax6,0:0.01:0.05);
colorbar(ax6); colormap(ax6, hot);

title(tile, sprintf("$r=%d$", r_opt), 'interpreter', 'latex');


%% PSDs of all the samples

FS_DHM = 115200;

psd_FOM = zeros(200, 2881);
psd_ROM = zeros(200, 2881);

overlap = round(0.75*split_size);

for i = 1:200
    [psd_seg, f_seg] = pwelch( ...
        EFOM(i,:), ...
        hann(split_size), overlap, split_size, FS_DHM);

    [psd_seg2, f_seg2] = pwelch( ...
        EROM{9}(i,:), ...
        hann(split_size), overlap, split_size, FS_DHM);

    % Normalize PSD
    psd_seg  = psd_seg ./ max(psd_seg);
    psd_seg2 = psd_seg2 ./ max(psd_seg2);

    psd_FOM(i,:)  = psd_seg;
    psd_ROM(i,:) = psd_seg2;
end
  
%%
figure;
for i = 1:200
    loglog(f_seg, psd_FOM(i,:), 'Color', [0.8 0.8 0.8], 'LineWidth', 0.5);
    hold on;
    loglog(f_seg, psd_ROM(i,:), 'Color', [0.8 0.8 0.8], 'LineWidth', 0.5);
end

hFOM = loglog(f_seg, mean(psd_FOM,1), 'r', 'LineWidth', 1);
hROM = loglog(f_seg, mean(psd_ROM,1), 'b', 'LineWidth', 1);

legend([hFOM, hROM], {'FOM mean PSD', 'ROM mean PSD'}, 'Location','best');
grid on;


%%
save(['/data/home/hyk049/Wave_Turbulence/' ...
    '0p10vpp_stoc_opinf_5datasets.mat'], 'EFOM', 'EROM');


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
