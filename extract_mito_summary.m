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

%% Whole-cell Vm and Cai (ci.txt)
% ci.txt columns (space-delimited, no header): time  vm  cai  ...
ci_data = readtable('ci.txt', 'FileType', 'text', 'Delimiter', ' ', 'ReadVariableNames', false);
t_ci   = ci_data{:, 1};
vm_ci  = ci_data{:, 2};
cai_ci = ci_data{:, 3};

figure('Name', 'Whole-Cell Vm and Cai', 'Position', [100 1100 900 500], 'Color', 'w');
subplot(2,1,1);
plot(t_ci, vm_ci, 'k', 'LineWidth', 1.2);
ylabel('V_m (mV)');
title('Whole-Cell Membrane Voltage and Ca_i');

subplot(2,1,2);
plot(t_ci, cai_ci, 'LineWidth', 1.2, 'Color', [0.00 0.45 0.74]);
ylabel('Ca_i (\muM)');
xlabel('Time (ms)');
%savefig('whole_cell_vm_cai_ATPD_0.1.fig')

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
%savefig('whole_cell_avg_psi_ca_ATPD_0.1.fig')

figure('Name', 'Whole-Cell ATP/ADP', 'Position', [100 100 900 950]);
yyaxis left;
plot(t, atp, 'LineWidth', 1.5);
ylabel('ATP (\muM)');
yyaxis right;
plot(t, adp, 'LineWidth', 1.5);
ylabel('ADP (\muM)');
%savefig('whole_cell_avg_ATP_ATPD_0.1.fig')

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

figure('Name', 'Prod Trace', 'Position', [1050 100 900 1050]);
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
%savefig('Prod_only_psi_ca_ATPD_0.1.fig')

figure('Name', 'Mito 0 Trace ATP/ADP', 'Position', [1050 100 900 1050]);
yyaxis left;
plot(t0, atp0, 'LineWidth', 1.5);
ylabel('ATP_0 (\muM)');
yyaxis right;
plot(t0, adp0, 'LineWidth', 1.5);
ylabel('ADP_0 (\muM)');

%savefig('Prod_only_ATP_ATPD_0.1.fig')
%Compare average ATP for simulation between producer and nonproduer crus
atp_grid = readmatrix('atp_full_grid.txt', 'FileType', 'text', 'Delimiter', '\t');
mask     = readmatrix('producer_mask.txt', 'FileType', 'text', 'Delimiter', '\t');

prod_cols    = find(mask == 1);
nonprod_cols = find(mask == 0);

atp_prod_avg    = mean(atp_grid(:, prod_cols), 2);
atp_nonprod_avg = mean(atp_grid(:, nonprod_cols), 2);
atp_diff        = atp_prod_avg - atp_nonprod_avg;

dt = 0.01;                                     % must match SAN_pace.cpp
output_interval = 100;                     
t_grid = (0:size(atp_grid,1)-1)' * output_interval * dt;   % ms

figure('Name', 'Producer vs Non-Producer ATP', 'Position', [1050 1300 900 400], 'Color', 'w');
plot(t_grid, atp_diff, 'LineWidth', 1.5, 'Color', [0.00 0.45 0.74]);
xlabel('Time (ms)');
ylabel('ATP_{producer} - ATP_{non-producer} (\muM)');
title('Producer vs. Non-Producer ATP Difference Over Time');
%savefig('ATP_Diff_ATPD_0.1.fig')
 