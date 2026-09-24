%% plot_mito_summary.m
% Quick sanity-check plot of mito_summary.txt
% Columns: time  avg_ATP  avg_ADP  avg_ca_mito  avg_psi_mito

data = readtable('mito_summary.txt', 'FileType', 'text', 'Delimiter', '\t');

t        = data.time;
atp      = data.avg_ATP;
adp      = data.avg_ADP;
ca_mito  = data.avg_ca_mito;
psi_mito = data.avg_psi_mito;

figure('Name', 'Mito Sanity Check', 'Position', [100 100 900 700]);

subplot(4,1,1);
plot(t, atp, 'LineWidth', 1.5);
ylabel('avg ATP (\muM)');
title('Mito Summary Averages Over Time (whole grid for ATP/ADP, producer CRUs for Ca/\Psi)');
grid on;

subplot(4,1,2);
plot(t, adp, 'LineWidth', 1.5, 'Color', [0.85 0.33 0.10]);
ylabel('avg ADP (\muM)');
grid on;

subplot(4,1,3);
plot(t, ca_mito, 'LineWidth', 1.5, 'Color', [0.47 0.67 0.19]);
ylabel('avg Ca_{mito} (\muM)');
grid on;

subplot(4,1,4);
plot(t, psi_mito, 'LineWidth', 1.5, 'Color', [0.49 0.18 0.56]);
ylabel('avg \Psi_{mito} (mV)');
xlabel('Time (ms)');
grid on;

%% Single-mito trace (mito0_trace.txt)
% Columns: time  cp0  J_uni0  jNCX_m0  ca_mito0  psi_mito0

data0 = readtable('mito0_trace.txt', 'FileType', 'text', 'Delimiter', '\t');

t0        = data0.time;
cp0       = data0.cp0;
J_uni0    = data0.J_uni0;
jNCX_m0   = data0.jNCX_m0;
ca_mito0  = data0.ca_mito0;
psi_mito0 = data0.psi_mito0;

figure('Name', 'Mito 0 Trace', 'Position', [1050 100 900 850]);

subplot(5,1,1);
plot(t0, cp0, 'LineWidth', 1.5);
ylabel('cp_0 (\muM)');
title('Mito 0 — Driving Ca and Fluxes Over Time');
grid on;

subplot(5,1,2);
plot(t0, J_uni0, 'LineWidth', 1.5, 'Color', [0.00 0.45 0.74]);
ylabel('J_{uni,0}');
grid on;

subplot(5,1,3);
plot(t0, jNCX_m0, 'LineWidth', 1.5, 'Color', [0.85 0.33 0.10]);
ylabel('J_{NCX,0}');
grid on;

subplot(5,1,4);
plot(t0, ca_mito0, 'LineWidth', 1.5, 'Color', [0.47 0.67 0.19]);
ylabel('Ca_{mito,0} (\muM)');
grid on;

subplot(5,1,5);
plot(t0, psi_mito0, 'LineWidth', 1.5, 'Color', [0.49 0.18 0.56]);
ylabel('\Psi_{mito,0} (mV)');
xlabel('Time (ms)');
grid on;

%% ATP linescan / kymograph (atp_linescan.txt)
% One row per output snapshot (every 100 steps), one column per CRU
% along x at fixed y = ny/2, z = nz/2 (mirrors Song's h_linescanBuffer).
% Tab-delimited, no header.

linescan = readmatrix('atp_linescan.txt', 'FileType', 'text', 'Delimiter', '\t');

figure('Name', 'ATP Linescan', 'Position', [100 850 900 400]);

imagesc(linescan');
set(gca, 'YDir', 'normal');
colorbar;
xlabel('Snapshot # (every 100 steps)');
ylabel('CRU position along x');
title('ATP_{cyto} Linescan (kymograph), y = ny/2, z = nz/2');