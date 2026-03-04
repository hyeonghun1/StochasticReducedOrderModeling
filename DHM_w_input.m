clear; clc;

%% 0.07vpp

input = '0p07vpp'; 

% Letters for the 10 files
% letters = ['a':'f', 'h':'k'];
letters = 'a':'j';

% letters = arrayfun(@(c) string(c), 'a':'i');
% letters(end+1) = "a_1127";

num_files = numel(letters);

base_path = '/disk/hyk049/DHM_new_1Dcenter/';

first_file = sprintf('%sQ_1D_%s_a.h5', base_path, input);
t = h5read(first_file, '/t')'; x = h5read(first_file, '/x')';
nx = length(x);

% Preallocate cell array for all Q_1D signals
Q_1D = cell(1, num_files);

% Read all files
for k = 1:num_files
    filename = base_path + "Q_1D_" + input + "_" + letters(k) + ".h5";
    Q_1D{k} = h5read(filename, '/Q_1D')';
end


%% Split samples

split_size = 5760;                 % # of snapshots per segment
num_segs = floor(length(t)/split_size);
tt = t(1:split_size);

num_signals = numel(Q_1D);

% Qstate_all: (nx, num_segs*num_signals, split_size)
Qstate_all = zeros(nx, num_segs*num_signals, split_size);

for s = 1:num_signals
    % Truncate each Q signal to an integer number of segments
    Q_all = Q_1D{s}(:, 1:split_size*num_segs);

    % Split each signal into segments and assign to Qstate_all
    for k = 1:num_segs
        idx_start = (k-1)*split_size + 1;
        idx_end   = k*split_size;

        % Compute column index range in Qstate_all
        col_idx = (s-1)*num_segs + k;

        % Assign the segment
        Qstate_all(:, col_idx, :) = Q_all(:, idx_start:idx_end);
    end
end

% Plot expected dynamics
[T, X] = meshgrid(tt, x);

figure;
surf(T, X, squeeze(mean(Qstate_all, 2)), 'EdgeColor', 'none');
xlabel('Time [s]'); ylabel('x [\mum]');
zlabel('Surface displacement [\mum]');
% xlim([0 0.06]); xticks(0:0.01:0.05); 
colorbar; colormap jet;


%% Filter signals

FS_DHM = 115200;
fc = 5e3;

% Design 4th-order Butterworth LPF
Wn = fc / (FS_DHM/2);
[b, a] = butter(4, Wn, 'low');

% Reshape to 2D for fast filtering
[nx, Nsamp, Nt] = size(Qstate_all);

Q_reshaped = reshape(Qstate_all, nx*Nsamp, Nt);

% Zero-phase filtering along time dimension
Q_filtered = filtfilt(b, a, Q_reshaped')';  % filtfilt works column-wise

% Reshape back to 3D
Qstate_filtered = reshape(Q_filtered, nx, Nsamp, Nt);


% Plot expected dynamics (filtered)
figure('color', 'w');
surf(T, X, squeeze(mean(Qstate_filtered, 2)), 'EdgeColor', 'none');
xlabel('Time [s]'); ylabel('x [\mum]');
zlabel('Surface displacement [\mum]');
colorbar; colormap jet;



%% FOM covariance
% CFOM = page_cov(Qstate_all, true);
% 
% %% plot empirical covariance
% 
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


%% Stochastic OpInf

% max r value
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


%% Filter signals
% FS_DHM = 115200;
% fc = 2e3;
% [high_freq, low_freq] = decompose_LPF(EFOM(100,:), FS_DHM, fc);
% 
% figure;
% plot(tt, low_freq)
%

%% PSD
% EFOM_traj = squeeze(mean(Qstate_all, 2));
% 
% window = hamming(1024);  % segment length
% noverlap = 512;           % 50% overlap
% nfft = 1024;              % FFT length
% 
% % [pxx, f] = pwelch(EFOM(100,:), window, noverlap, nfft, FS_DHM);
% [pxx2, f2] = pwelch(EFOM_traj(100,:), window, noverlap, nfft, FS_DHM);
% 
% figure; hold on;
% % plot(f, log10(pxx));
% plot(f2, log10(pxx2));
% xlim([10, fc+3000]);
% xlabel('Frequency (Hz)');
% ylabel('PSD (dB/Hz)');


%%
% Test with different number of noise samples

% u_train = (BC_1 + BC_2)/2 ;

f_input = 7e6;
u_amp = 1;
u_train = u_amp*cos(2*pi*f_input*tt);
% figure;
% plot(tt, u_train)


%% Test with fixed u_amp, sigma + regularization

A_reg = logspace(-1, 3, 10);
B_reg = logspace(0, 8, 10);
N_reg = logspace(0, 5, 10);

r_opt = 8;
Vr_temp = V(:, 1:r_opt);      % ROM basis
xr0 = Vr_temp' * EFOM(:,1);   % ROM initial consdtion

Q_train_temp = pagemtimes(Vr_temp', Qstate_all);  % R^{r x L x s}
E_train = squeeze(mean(Q_train_temp, 2));         % R^{r x s}
C_train = page_cov(Q_train_temp, true);           % R^{r x r x s}
    
isbilinear = true;

if isbilinear
    outerLoop = 1:length(N_reg);
    ROM_form = 'non-autonomous bilinear ROM';
else
    outerLoop = 1:length(A_reg); 
    ROM_form = 'non-autonomous linear ROM';
end

% Initialization
EROM = cell(1, length(A_reg));
CROM = cell(1, length(B_reg));
E_error = zeros(length(A_reg), length(B_reg));
C_error = zeros(length(A_reg), length(B_reg));


rng(seed_test)
for ii = outerLoop
    tic
    disp(ii + "th A_reg")

    EROM{ii} = cell(1, length(B_reg));
    CROM{ii} = cell(1, length(B_reg));    
    
    % rng(seed_test)

    for jj = 1:length(B_reg)
    
    % rng(seed_test)

    if isbilinear
        lambda = [1, B_reg(jj), N_reg(ii)];
    else
        lambda = [A_reg(ii), B_reg(jj), 1];
    end

    % Drift OpInf
    [Ehatr, Ahatr, Bhatr, Nhatr] = infer_drift_u(E_train, u_train, h, isbilinear, ...
                                                    lambda);
    % Diffusion OpInf
    [Mhatr, Khatr] = infer_diffusion_u(C_train, u_train, h, Ahatr, Nhatr);
    
    sigma = 20;
    % rng(seed_test)
    % ROM function using Wiener process
    fhatr = @(x0, u, L) ...
            (Ehatr - h*Ahatr - h*Nhatr*u) \ ...
            (x0 + h*Bhatr*u + sqrt(h)*Mhatr*sigma*randn(size(Mhatr,2), L)) ;
    
    % ROM simulation w/ testing I.C. across L samples
    % Eopinf: R^{r x s}, Copinf: R^{r x r x s} 
    [Eopinf, Copinf] = compute_model_u(fhatr, Vr_temp, xr0, u_train, L, L);

    % Reconstruct FOM dimension
    E_recon = Vr_temp * Eopinf;  % R^{n x s}
     
    EROM{ii}{jj} = E_recon;
    CROM{ii}{jj} = Copinf;
    
    E_error(ii,jj) = norm(E_train - Eopinf, 'fro') / norm(E_train, 'fro');
    C_error(ii,jj) = page_norm(C_train - Copinf) / page_norm(C_train);
    end
    
    fprintf('Elapsed time: %.2f seconds.\n', toc);
end

error.E_error = E_error;
error.C_error = C_error;


%%
save(['/disk/hyk049/WT_RomFit/0p07/' ...
    '0p07vpp_r=8_ABN_amp1_sigma20.mat'], 'EFOM', 'EROM_opt');


%% Plot fitted result

[TT, X] = meshgrid(tt, x);
[i_min, j_min] = find(E_error == min(E_error(:)));
EROM_opt = EROM{i_min}{j_min};

rel_error = abs(EFOM - EROM_opt) ./ abs(EFOM);

figure('Color','w');
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

title(tile, sprintf("$r=%d$", r_opt), 'Interpreter','latex');

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


%% Test with different input amp, noise amplitudes

A_reg = logspace(-1, 3, 10);
B_reg = logspace(0, 8, 10);
N_reg = logspace(0, 5, 10);

r_opt = 18;
Vr_temp = V(:, 1:r_opt);

xr0 = Vr_temp' * EFOM(:,1);  % ROM initial condition

Q_train_temp = pagemtimes(Vr_temp', Qstate_all);  % R^{r x L x s}
E_train = squeeze(mean(Q_train_temp, 2));         % R^{r x s}
C_train = page_cov(Q_train_temp, true);           % R^{r x r x s}
    
isbilinear = false;

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

EROM1  = EROM_all.amp100.sigma1;
EROM10 = EROM_all.amp100.sigma10;
EROM15 = EROM_all.amp100.sigma15;
EROM20 = EROM_all.amp100.sigma20;

rel1  = abs(EFOM - EROM1)  ./ abs(EFOM);
rel10 = abs(EFOM - EROM10) ./ abs(EFOM);
rel15 = abs(EFOM - EROM15) ./ abs(EFOM);
rel20 = abs(EFOM - EROM20) ./ abs(EFOM);

figure('Color','w');
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

function [stationary, drift] = decompose_LPF(x, fs, fc)
    % 3rd-order Butterworth low-pass filter
    [b, a] = butter(3, fc/(fs/2), 'low');
    
    % Zero-phase filtering (filtfilt)
    drift = filtfilt(b, a, x);
    
    % Stationary = original - drift
    stationary = x - drift;
end