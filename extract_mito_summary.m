%% plot_mito_summary.m
% Quick sanity-check plot of mito_summary.txt
% Columns: time  avg_ATP  avg_ADP  avg_ca_mito  avg_psi_mito  avg_cai  avg_cp_prod
clear all; close all;

data = readtable('mito_summary.txt', 'FileType', 'text', 'Delimiter', '\t');

t         = data.time;
atp       = data.avg_ATP;
adp       = data.avg_ADP;
ca_mito   = data.avg_ca_mito;
psi_mito  = data.avg_psi_mito;
avg_cai   = data.avg_cai;       % whole-cell average cytosolic Ca
avg_cp    = data.avg_cp_prod;   % cleft Ca, averaged over producer CRUs only (mito Ca input)

figure('Name', 'Whole-Cell Averages Ca and Mito', 'Position', [100 100 900 950]);

subplot(4,1,1);
plot(t, avg_cai, 'LineWidth', 1.5, 'Color', [0.00 0.45 0.74]);
ylabel('avg Ca_i (\muM)');

subplot(4,1,2);
plot(t, avg_cp, 'LineWidth', 1.5, 'Color', [0.93 0.69 0.13]);
ylabel('avg Ca_{cleft,prod} (\muM)');


subplot(4,1,3);
plot(t, ca_mito, 'LineWidth', 1.5, 'Color', [0.47 0.67 0.19]);
ylabel('avg Ca_{mito} (\muM)');


subplot(4,1,4);
plot(t, psi_mito, 'LineWidth', 1.5, 'Color', [0.49 0.18 0.56]);
ylabel('avg \Psi_{mito} (mV)');
xlabel('Time (ms)');

figure('Name', 'Whole-Cell ATP/ADP', 'Position', [100 100 900 950]);
yyaxis left;
plot(t, atp, 'LineWidth', 1.5);
ylabel('ATP (\muM)');
yyaxis right;
plot(t, adp, 'LineWidth', 1.5);
ylabel('ADP (\muM)');

%% Single-mito trace (mito0_trace.txt)
% Columns: time  cp0  cai0  J_uni0  jNCX_m0  ca_mito0  psi_mito0  atp0  adp0  prod0  consum0  diff0

data0 = readtable('mito0_trace.txt', 'FileType', 'text', 'Delimiter', '\t');

t0        = data0.time;
cp0       = data0.cp0;
cai0      = data0.cai0;
J_uni0    = data0.J_uni0;
jNCX_m0   = data0.jNCX_m0;
ca_mito0  = data0.ca_mito0;
psi_mito0 = data0.psi_mito0;
atp0      = data0.atp0;
adp0      = data0.adp0;
prod0     = data0.prod0;    % local ATP production rate (VATPase) at mito0's producer CRU
consum0   = data0.consum0;  % local ATP consumption rate (VATP_consum) at that same CRU
diff0     = data0.diff0;    % net diffusive term (J_ATP_D) at that same CRU

figure('Name', 'Mito 0 Trace', 'Position', [1050 100 900 1050]);

subplot(4,1,1);
plot(t0, cp0, 'LineWidth', 1.5);
ylabel('cp_0 (\muM)');
title('Mito 0 — Driving Ca, Fluxes, and ATP/ADP Over Time');

subplot(4,1,2);
plot(t0, cai0, 'LineWidth', 1.5, 'Color', [0.00 0.45 0.74]);
ylabel('cai_0 (\muM)');

% subplot(7,1,3);
% plot(t0, J_uni0, 'LineWidth', 1.5, 'Color', [0.00 0.45 0.74]);
% hold on;
% plot(t0, jNCX_m0, 'LineWidth', 1.5, 'Color', [0.85 0.33 0.10]);
% hold off;
% ylabel('J_{uni,0} / J_{NCX,0}');
% legend({'J_{uni,0}', 'J_{NCX,0}'}, 'Location', 'best');

subplot(4,1,3);
plot(t0, psi_mito0, 'LineWidth', 1.5, 'Color', [0.49 0.18 0.56]);
ylabel('\Psi_{mito,0} (mV)');

subplot(4,1,4);
plot(t0, ca_mito0, 'LineWidth', 1.5, 'Color', [0.47 0.67 0.19]);
ylabel('Ca_{mito,0} (\muM)');

figure('Name', 'Mito 0 Trace ATP/ADP', 'Position', [1050 100 900 1050]);
yyaxis left;
plot(t0, atp0, 'LineWidth', 1.5);
ylabel('ATP_0 (\muM)');
yyaxis right;
plot(t0, adp0, 'LineWidth', 1.5);
ylabel('ADP_0 (\muM)');

% subplot(7,1,7);
% plot(t0, prod0, 'LineWidth', 1.5, 'Color', [0.00 0.60 0.20]);
% hold on;
% plot(t0, consum0, 'LineWidth', 1.5, 'Color', [0.80 0.10 0.10]);
% plot(t0, diff0, 'LineWidth', 1.5, 'Color', [0.30 0.30 0.30]);
% hold off;
% ylabel('rate (\muM/ms)');
% xlabel('Time (ms)');
% legend({'production (V_{ATPase})', 'consumption (V_{ATP,consum})', 'net diffusion (J_{ATP,D})'}, 'Location', 'best');
% title('Competing ATP rate terms at mito 0''s producer CRU');


%% ATP linescan (atp_linescan.txt)
% One row per output snapshot (every 100 steps), one column per CRU
% along y at fixed x = nx/2, z = nz/2. Tab-delimited, no header.

 
%data = readmatrix('atp_linescan.txt', 'FileType', 'text', 'Delimiter', '\t');
 
%figure('Position', [100 850 900 400], 'Color', 'w');
 
%[rows, cols] = size(data);
%scale_factor = 3;   % try 2-5; higher = smoother
%[X, Y]   = meshgrid(1:cols, 1:rows);
%[Xq, Yq] = meshgrid(linspace(1, cols, cols*scale_factor), linspace(1, rows, rows*scale_factor));
%data_smooth = interp2(X, Y, data, Xq, Yq, 'spline');
 
%imagesc(data_smooth');
%colormap('jet');
% caxis([4955, 5000]);   % set once you know a good comparison range across runs; auto-scaled for now
%set(gca, 'YTick', [], 'XTick', [], 'XTickLabel', [], 'Color', 'w', 'XColor', 'k', 'YColor', 'k');
%gca_ax = gca; gca_ax.XAxis.Visible = 'off';
%ylabel(sprintf('%.1f \\mum', rows * 0.9));   % ny * l_T_atp, physical extent along y
%cb = colorbar;
%cb.Label.String = 'ATP_{cyto} (\muM)';
%cb.Color = 'k';
%box off;
%% ATP linescan (atp_linescan.txt), styled to match Haibo's Ca linescan script
% One row per output snapshot (every 100 steps), one column per CRU
% along y at fixed x = nx/2, z = nz/2. Tab-delimited, no header.
% data layout matches Haibo's: rows = time, columns = space.

data = readmatrix('atp_linescan.txt', 'FileType', 'text', 'Delimiter', '\t');

figure('Position', [100 850 900 400], 'Color', 'w');

[rows, cols] = size(data);
scale_factor = 3;   % try 2-5; higher = smoother
[X, Y]   = meshgrid(1:cols, 1:rows);
[Xq, Yq] = meshgrid(linspace(1, cols, cols*scale_factor), linspace(1, rows, rows*scale_factor));
data_smooth = interp2(X, Y, data, Xq, Yq, 'spline');

imagesc(data_smooth');
colormap('jet');
% caxis([4955, 5000]);   % set once you know a good comparison range across runs; auto-scaled for now
set(gca, 'YTick', [], 'XTick', [], 'XTickLabel', [], 'Color', 'w', 'XColor', 'k', 'YColor', 'k');
gca_ax = gca; gca_ax.XAxis.Visible = 'off';
ylabel(sprintf('%.1f \\mum', rows * 0.9));   % ny * l_T_atp, physical extent along y
cb = colorbar;
cb.Label.String = 'ATP_{cyto} (\muM)';
cb.Color = 'k';
box off;

%% Same linescan, but as deviation from each snapshot's own spatial mean
% Raw ATP_cyto sits on a slowly drifting ~5000 uM baseline that swamps the
% actual spatial signal (~0.01-0.3 uM) visually. Subtracting each row's
% (snapshot's) own mean removes that common drift and rescales the color
% axis to just the spatial structure -- same trick as dF/F0 normalization
% in Ca imaging.

data_dev = data - mean(data, 2);   % subtract each snapshot's row-mean

[rows_d, cols_d] = size(data_dev);
[Xd, Yd]   = meshgrid(1:cols_d, 1:rows_d);
[Xqd, Yqd] = meshgrid(linspace(1, cols_d, cols_d*scale_factor), linspace(1, rows_d, rows_d*scale_factor));
data_dev_smooth = interp2(Xd, Yd, data_dev, Xqd, Yqd, 'spline');

figure('Position', [100 1300 900 400], 'Color', 'w');

imagesc(data_dev_smooth');
colormap('jet');
clim_val = max(abs(data_dev(:)));
caxis([-clim_val, clim_val]);   % symmetric around zero so the diverging colormap is centered on "no deviation"
set(gca, 'YTick', [], 'XTick', [], 'XTickLabel', [], 'Color', 'w', 'XColor', 'k', 'YColor', 'k');
gca_ax = gca; gca_ax.XAxis.Visible = 'off';
ylabel(sprintf('%.1f \\mum', rows_d * 0.9));
cb = colorbar;
cb.Label.String = 'ATP_{cyto} deviation from mean (\muM)';
cb.Color = 'k';
box off;
title('Same linescan, deviation from each snapshot''s spatial mean', 'Color', 'k');



 