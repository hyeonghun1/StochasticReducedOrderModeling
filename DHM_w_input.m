clear; clc;

%% 0.25vpp

t = h5read('/disk/hyk049/DHM_new_1Dcenter/Q_1D_0p25vpp_a.h5', '/t')';
x = h5read('/disk/hyk049/DHM_new_1Dcenter/Q_1D_0p25vpp_a.h5', '/x')';
nx = length(x);

Q_a_1D = h5read('/disk/hyk049/DHM_new_1Dcenter/Q_1D_0p25vpp_a.h5', '/Q_1D')';
Q_b_1D = h5read('/disk/hyk049/DHM_new_1Dcenter/Q_1D_0p25vpp_b.h5', '/Q_1D')';
Q_c_1D = h5read('/disk/hyk049/DHM_new_1Dcenter/Q_1D_0p25vpp_c.h5', '/Q_1D')';
Q_d_1D = h5read('/disk/hyk049/DHM_new_1Dcenter/Q_1D_0p25vpp_d.h5', '/Q_1D')';
Q_e_1D = h5read('/disk/hyk049/DHM_new_1Dcenter/Q_1D_0p25vpp_e.h5', '/Q_1D')';
Q_f_1D = h5read('/disk/hyk049/DHM_new_1Dcenter/Q_1D_0p25vpp_f.h5', '/Q_1D')';
Q_g_1D = h5read('/disk/hyk049/DHM_new_1Dcenter/Q_1D_0p25vpp_g.h5', '/Q_1D')';
Q_h_1D = h5read('/disk/hyk049/DHM_new_1Dcenter/Q_1D_0p25vpp_h.h5', '/Q_1D')';
Q_i_1D = h5read('/disk/hyk049/DHM_new_1Dcenter/Q_1D_0p25vpp_i.h5', '/Q_1D')';
Q_j_1D = h5read('/disk/hyk049/DHM_new_1Dcenter/Q_1D_0p25vpp_j.h5', '/Q_1D')';


%% 10 datasets centerline case

split_size = 5760; % # of snapshots with
% one segment of temporal independence

num_segs = floor(length(t)/split_size);
tt = t(1:split_size);

Q_split_a = cell(1, num_segs);
Q_split_b = cell(1, num_segs);
Q_split_c = cell(1, num_segs);
Q_split_d = cell(1, num_segs);
Q_split_e = cell(1, num_segs);
Q_split_f = cell(1, num_segs);
Q_split_g = cell(1, num_segs);
Q_split_h = cell(1, num_segs);
Q_split_i = cell(1, num_segs);
Q_split_j = cell(1, num_segs);

Q_all_a = Q_a_1D(:,1:split_size*num_segs);
Q_all_b = Q_b_1D(:,1:split_size*num_segs);
Q_all_c = Q_c_1D(:,1:split_size*num_segs);
Q_all_d = Q_d_1D(:,1:split_size*num_segs);
Q_all_e = Q_e_1D(:,1:split_size*num_segs);
Q_all_f = Q_f_1D(:,1:split_size*num_segs);
Q_all_g = Q_g_1D(:,1:split_size*num_segs);
Q_all_h = Q_h_1D(:,1:split_size*num_segs);
Q_all_i = Q_i_1D(:,1:split_size*num_segs);
Q_all_j = Q_j_1D(:,1:split_size*num_segs);


for k=1:num_segs
    idx_start = (k-1)*split_size + 1;
    idx_end = k*split_size;

    Q_split_a{k} = Q_all_a(:, idx_start:idx_end);
    Q_split_b{k} = Q_all_b(:, idx_start:idx_end);
    Q_split_c{k} = Q_all_c(:, idx_start:idx_end);
    Q_split_d{k} = Q_all_d(:, idx_start:idx_end);
    Q_split_e{k} = Q_all_e(:, idx_start:idx_end);
    Q_split_f{k} = Q_all_f(:, idx_start:idx_end);
    Q_split_g{k} = Q_all_g(:, idx_start:idx_end);
    Q_split_h{k} = Q_all_h(:, idx_start:idx_end);
    Q_split_i{k} = Q_all_i(:, idx_start:idx_end);
    Q_split_j{k} = Q_all_j(:, idx_start:idx_end);
end

Qstate_a = zeros(nx, num_segs, split_size);
Qstate_b = zeros(nx, num_segs, split_size);
Qstate_c = zeros(nx, num_segs, split_size);
Qstate_d = zeros(nx, num_segs, split_size);
Qstate_e = zeros(nx, num_segs, split_size);
Qstate_f = zeros(nx, num_segs, split_size);
Qstate_g = zeros(nx, num_segs, split_size);
Qstate_h = zeros(nx, num_segs, split_size);
Qstate_i = zeros(nx, num_segs, split_size);
Qstate_j = zeros(nx, num_segs, split_size);

for i=1:num_segs
    Qstate_a(:,i,:) = Q_split_a{i};
    Qstate_b(:,i,:) = Q_split_b{i};
    Qstate_c(:,i,:) = Q_split_c{i};
    Qstate_d(:,i,:) = Q_split_d{i};
    Qstate_e(:,i,:) = Q_split_e{i};
    Qstate_f(:,i,:) = Q_split_f{i};
    Qstate_g(:,i,:) = Q_split_g{i};
    Qstate_h(:,i,:) = Q_split_h{i};
    Qstate_i(:,i,:) = Q_split_i{i};
    Qstate_j(:,i,:) = Q_split_j{i};
end


%%
Qstate_all = cat(2, Qstate_a, Qstate_b, Qstate_c, Qstate_d, Qstate_e, ...
                Qstate_f, Qstate_g, Qstate_h, Qstate_i, Qstate_j);

[T, X] = meshgrid(tt, x);

figure;
surf(T, X, squeeze(mean(Qstate_all, 2)), 'EdgeColor', 'none');
xlabel('Time [s]'); ylabel('x [\mum]');
zlabel('Surface displacement [\mum]');
% xlim([0 0.06]); xticks(0:0.01:0.05); 
colorbar; colormap jet;

%%
Q_subspace = squeeze(mean(Qstate_all, 2));
min(Q_subspace(:))


%%
[T, X] = meshgrid(tt, x);

Z = cell(1,10);
Z{1} = squeeze(mean(Qstate_a, 2));
Z{2} = squeeze(mean(Qstate_b, 2));
Z{3} = squeeze(mean(Qstate_c, 2));
Z{4} = squeeze(mean(Qstate_d, 2));
Z{5} = squeeze(mean(Qstate_e, 2));
Z{6} = squeeze(mean(Qstate_f, 2));
Z{7} = squeeze(mean(Qstate_g, 2));
Z{8} = squeeze(mean(Qstate_h, 2));
Z{9} = squeeze(mean(Qstate_i, 2));
Z{10} = squeeze(mean(Qstate_j, 2));

% Compute common color limits (exclude last if not plotting it)
zmin = min([Z{1}(:); Z{2}(:); Z{3}(:); Z{4}(:); Z{5}(:); ...
            Z{6}(:); Z{7}(:); Z{8}(:); Z{9}(:); Z{10}(:)]);
zmax = max([Z{1}(:); Z{2}(:); Z{3}(:); Z{4}(:); Z{5}(:); ...
            Z{6}(:); Z{7}(:); Z{8}(:); Z{9}(:); Z{10}(:)]);
figure;
tiledlayout(2, 5,'Padding','compact','TileSpacing','compact');

for k = 1:10
    ax = nexttile;
    surf(ax, T, X, Z{k}, 'EdgeColor', 'none');
    xlabel(ax,'Time [s]'); ylabel(ax,'x [\mum]');
    zlabel(ax,'Surface displacement [\mum]');
    xlim(ax,[0 0.06]); xticks(ax,0:0.01:0.05);
    caxis(ax,[zmin zmax]);
    view(ax,3);
    colormap(ax, jet);
    colorbar(ax);
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


%% Get POD basis

disp("Computing subspace Vr...")
[V,S,~] = svd(reshape(Qstate_all, nx, []), "econ");
disp("Subspace computation complete!")


%% plot singular values &
% Projection errors

normalized_sigma = diag(S)/max(diag(S));

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


%% plot projections V*V^T*Q
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
% xlim([0 0.06]); xticks(0:0.01:0.05);
view(2)   % optional: top view (remove if you want 3D)

% Plot ROMs
for k = 1:numel(r_list)
    r = r_list(k);

    Vrtemp = V(:, 1:r);
    Q_proj = Vrtemp * (Vrtemp' * reshape(Qstate_all, nx, []));
    Q_proj = reshape(Q_proj, nx, size(Qstate_all,2), size(Qstate_all,3));

    QROM_mean = squeeze(mean(Q_proj, 2));

    ax(k+1) = nexttile;
    s(k+1) = surf(T, X, QROM_mean, 'EdgeColor','none');
    title(['Projection (r = ' num2str(r) ')'])
    xlabel('Time [s]'); ylabel('x [\mum]');
    % xlim([0 0.06]); xticks(0:0.01:0.05);
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


%% Stochastic OpInf

% max r value
rmax = 20;

h = double(tt(2) - tt(1));
[~, L, s] = size(Qstate_all);
isbilinear = false;

seed_test = 42;

% FOM expectation (average) over all samples
EFOM = squeeze(mean(Qstate_all, 2)); 

% Testing noise sample numbers
L_test = L * [1];
% L_test = L * [1, 10, 100, 1000];


%%
BC_1 = EFOM(1,:);
BC_2 = EFOM(end,:);
u_train = (BC_1 + BC_2)/2 ;

figure;
plot(BC_1, 'DisplayName', '1'); 
hold on
plot(BC_2, 'DisplayName', 'end');
plot(u_train, 'DisplayName', 'u_{input}');

legend show
grid on

%%
% Test with different number of noise samples

% u_train = (BC_1 + BC_2)/2 ;

f_input = 7e6;
u_amp = 10;
u_train = u_amp*cos(2*pi*f_input*tt);
figure;
plot(tt, u_train)

%%

A_reg = logspace(-1, 2, 10);
B_reg = logspace(2, 6, 10);
N_reg = logspace(0, 3, 10);

EROM = cell(1, length(A_reg));
CROM = cell(1, length(B_reg));
E_error = zeros(length(A_reg), length(B_reg));
C_error = zeros(length(A_reg), length(B_reg));

r_opt = 17;
Vr_temp = V(:, 1:r_opt);

% Project the FOM observations
Q_train_temp = pagemtimes(Vr_temp', Qstate_all);  % R^{r x L x s}

% Estimate mean/covariance of the reduced observations
% mean of reduced states across L realizations
E_train = reshape(squeeze(mean(Q_train_temp, 2)), r_opt, s); % R^{r x s}

% Covariance of reduced states
C_train = page_cov(Q_train_temp, true); % R^{r x r x s}
    
isbilinear = false;

for ii = 1:length(A_reg)
    tic
    disp(ii + "th A_reg")

    EROM{ii} = cell(1, length(B_reg));
    CROM{ii} = cell(1, length(B_reg));    

    for jj = 1:length(B_reg)
    
    rng(seed_test)
    lambda = [A_reg(ii), B_reg(jj), 1];
    
    % Drift OpInf
    [Ehatr, Ahatr, Bhatr, Nhatr] = infer_drift_u(E_train, u_train, h, isbilinear, lambda);
    
    % Diffusion OpInf
    [Mhatr, Khatr] = infer_diffusion_u(C_train, u_train, h, Ahatr, Nhatr);
    
    sigma = 10;
    % ROM function using Wiener process
    fhatr = @(x0, u, L) ...
            (Ehatr - h*Ahatr - h*Nhatr*u) \ ...
            (x0 + h*Bhatr*u + sqrt(h)*Mhatr*sigma*randn(size(Mhatr,2), L)) ;
    
    % ROM initial consdtion
    xr0 = Vr_temp' * EFOM(:,1);
    
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
    
    toc
end

error.E_error = E_error;
error.C_error = C_error;



%%
save(['/disk/hyk049/WT_RomFit/0p25/' ...
    '0p25vpp_stoc_opinf_10datasets_B_10sigma.mat'], 'EFOM', 'EROM_opt');


%%
[TT, X] = meshgrid(tt, x);

EROM_opt = EROM{9}{5};

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



%%
% Test with different number of noise samples

u_train = (BC_1 + BC_2)/2 ;

for j=1:length(L_test)

    EROM{j} = cell(1,rmax); CROM{j} = cell(1,rmax);    
    E_error{j} = zeros(1,rmax); C_error{j} = zeros(1,rmax);

    tic
    for ii = 1:rmax
        disp("ROM dimension " + ii)
        rng(seed_test)
    
        Vr_temp = V(:, 1:ii);
    
        % Project the FOM observations
        Q_train_temp = pagemtimes(Vr_temp', Qstate_all);  % R^{r x L x s}
    
        % Estimate mean/covariance of the reduced observations
        % mean of reduced states across L realizations
        E_train = reshape(squeeze(mean(Q_train_temp, 2)), ii, s); % R^{r x s}
    
        % Covariance of reduced states
        C_train = page_cov(Q_train_temp, true); % R^{r x r x s}
    
        isbilinear = false;
        lambda = [1, 1e3];
        % Drift OpInf
        [Ehatr, Ahatr, Bhatr, Nhatr] = infer_drift_u(E_train, u_train, h, isbilinear, lambda);
    
        % Diffusion OpInf
        [Mhatr, Khatr] = infer_diffusion_u(C_train, u_train, h, Ahatr, Nhatr);
    
        sigma = 1;
        % ROM function using Wiener process
        fhatr = @(x0, u, L) ...
            (Ehatr - h*Ahatr - h*Nhatr*u) \ ...
            (x0 + h*Bhatr*u + sqrt(h)*Mhatr*sigma*randn(size(Mhatr,2), L)) ;

        % fhatr = @(x0, L) ...
        %     (Ehatr - h*Ahatr) \ ...
        %     (x0 + sqrt(h)*Mhatr*sigma*randn(size(Mhatr,2), L)) ;
    
        % ROM initial consdtion
        xr0 = Vr_temp' * EFOM(:,1);
    
        % ROM simulation w/ testing I.C.
        % This simulates ROM across L samples
        % Eopinf: R^{r x s}, Copinf: R^{r x r x s} 
        [Eopinf, Copinf] = compute_model_u(fhatr, Vr_temp, ...
                                            xr0, u_train, L, L_test(j));
        % Reconstruct FOM dimension
        E_recon = Vr_temp * Eopinf;  % R^{n x s}
    
        % Vr*C*Vr^T (R^{n x n x s})
        % C_recon = pagemtimes( pagemtimes(Vr_temp, Copinf), Vr_temp' );
     
        EROM{j}{ii} = E_recon;
        CROM{j}{ii} = Copinf;
    
        E_error{j}(ii) = norm(E_train - Eopinf, 'fro') / norm(E_train, 'fro');
        C_error{j}(ii) = page_norm(C_train - Copinf) / page_norm(C_train);
    end
    toc

end

error.E_error = E_error;
error.C_error = C_error;

%% plot errors
r = linspace(1, rmax, rmax);

figure('Color', 'w');

for j=1:length(L_test)
% expectation
subplot(1,2,1)
plot(r, error.E_error{j}, '-o', 'linewidth', 1.5, ...
    'MarkerFaceColor','auto'); hold on;
set(gca, 'YScale', 'log')
xlabel("Reduced dimension $r$", 'interpreter', 'latex')
ylabel("$e_E$", 'interpreter', 'latex')
title("Expectation error")
% axis([1 rmax min(error.E_error)*0.9 max(error.E_error)*1.1])
grid on

% Covariance
subplot(1,2,2)
plot(r, error.C_error{j}, '-o', 'LineWidth', 1.5, ...
    'MarkerFaceColor','auto'); hold on;
set(gca, 'YScale', 'log')
xlabel("Reduced dimension $r$", 'interpreter', 'latex')
ylabel("$e_C$", 'interpreter', 'latex')
title("Covariance error")
% axis([1 rmax min(error.C_error)*0.9 max(error.C_error)*1.1])
grid on

legend(sprintf("L=%d", L_test(1)), sprintf("L=%d", L_test(2)));

end


%%
[TT, X] = meshgrid(tt, x);

r_opt = 13;
EROM_opt = EROM{1}{r_opt};

figure('Color','w');
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
title(ax2,"Expected dynamics of" + newline + "stochastic OpInf ROM fitting" ...
      + newline + "(L=170,000)");
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


%%
LL = 170000;  % testing # of noise samples
rr = 18;
Vr_temp = V(:, 1:rr);

h = double(tt(2) - tt(1));
[~, L, s] = size(Qstate_all);
isbilinear = false;

% FOM expectation (average) over all samples
EFOM = squeeze(mean(Qstate_all, 2)); 

seed_test = 42;
rng(seed_test);

% Project the FOM observations
Q_train_temp = pagemtimes(Vr_temp', Qstate_all);  % R^{r x L x s}

% Estimate mean/covariance of the reduced observations
% mean of reduced states across L realizations
E_train = reshape(squeeze(mean(Q_train_temp, 2)), rr, s); % R^{r x s}

% Covariance of reduced states
C_train = page_cov(Q_train_temp, true); % R^{r x r}

% Drift OpInf
[Ehatr, Ahatr, Nhatr] = infer_drift(E_train, h, isbilinear, s);

% Diffusion OpInf
[Mhatr, Khatr] = infer_diffusion(C_train, h, Ahatr, Nhatr);

sigma = 5;
% ROM function using Wiener process
fhatr = @(x0, L) (Ehatr - h*Ahatr) \ ...
                    (x0 + sqrt(h)*Mhatr*sigma*randn(size(Mhatr,2), L)) ;

% ROM initial consdtion
xr0 = Vr_temp' * EFOM(:,1);

% ROM simulation w/ testing I.C.
% This simulates ROM across L samples
% Eopinf: R^{r x s}, Copinf: R^{r x r x s}
[Eopinf, Copinf, f1opinf, f2opinf] = compute_model(fhatr, Vr_temp, xr0, ...
                                                    s, L, LL);
% Reconstruct FOM dimension
E_recon = Vr_temp * Eopinf;  % R^{n x s}

% Vr*C*Vr^T (R^{n x n x s})
C_recon = pagemtimes( pagemtimes(Vr_temp, Copinf), Vr_temp' );


%%
save(['/data/home/hyk049/Wave_Turbulence/' ...
    '0p15vpp_stoc_opinf_10datasets_input.mat'], 'EFOM', 'EROM_opt');


%%
EROM_opt = E_recon;

[TT, X] = meshgrid(tt, x);

figure('Color','w');
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
title(ax2,"Expected dynamics of" + newline + "stochastic OpInf ROM fitting" ...
      + newline + "(L=170,000)");
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

title(tile, sprintf("$r=%d$", rr), 'interpreter', 'latex');

zmin = min([s1.ZData(:); s2.ZData(:); s4.ZData(:); s5.ZData(:)]);
zmax = max([s1.ZData(:); s2.ZData(:); s4.ZData(:); s5.ZData(:)]);
set([ax1 ax2 ax4 ax5], 'CLim', [zmin zmax])

zlim(ax1, [zmin zmax])
zlim(ax2, [zmin zmax])
zlim(ax4, [zmin zmax])
zlim(ax5, [zmin zmax])



%% Functions
function y = page_norm(A)
    y = sum(abs(A).^2, 'all');
    y = sqrt(y);
end