%% plot_mito_summary.m
% Quick sanity-check plot of mito_summary.txt
% Columns: time  avg_ATP  avg_ADP  avg_ca_mito  avg_psi_mito  avg_cai  avg_cp_prod
clear all; close all;
run_label = '4_4_spont_ikatp_0';  % Label for graphs + biomarkers file
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
savefig(['whole_cell_vm_cai_', run_label, '.fig'])

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
savefig(['whole_cell_avg_psi_ca_', run_label, '.fig'])

figure('Name', 'Whole-Cell ATP/ADP', 'Position', [100 100 900 950]);
yyaxis left;
plot(t, atp, 'LineWidth', 1.5);
ylabel('ATP (\muM)');
yyaxis right;
plot(t, adp, 'LineWidth', 1.5);
ylabel('ADP (\muM)');
savefig(['whole_cell_avg_ATP_ATPD', run_label, '.fig'])

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
% 
% % subplot(7,1,3);
% % plot(t0, J_uni0, 'LineWidth', 1.5, 'Color', [0.00 0.45 0.74]);
% % hold on;
% % plot(t0, jNCX_m0, 'LineWidth', 1.5, 'Color', [0.85 0.33 0.10]);
% % hold off;
% % ylabel('J_{uni,0} / J_{NCX,0}');
% % legend({'J_{uni,0}', 'J_{NCX,0}'}, 'Location', 'best');
% 
subplot(4,1,3);
plot(t0, psi_mito0, 'LineWidth', 1.5, 'Color', [0.49 0.18 0.56]);
ylabel('\Psi_{mito,0} (mV)');

subplot(4,1,4);
plot(t0, ca_mito0, 'LineWidth', 1.5, 'Color', [0.47 0.67 0.19]);
ylabel('Ca_{mito,0} (\muM)');
savefig(['Prod_only_psi_ca_ATPD_', run_label, '.fig'])
% 
% figure('Name', 'Mito 0 Trace ATP/ADP', 'Position', [1050 100 900 1050]);
% yyaxis left;
% plot(t0, atp0, 'LineWidth', 1.5);
% ylabel('ATP_0 (\muM)');
% yyaxis right;
% plot(t0, adp0, 'LineWidth', 1.5);
% ylabel('ADP_0 (\muM)');

% % %savefig(['Prod_only_ATP_ATPD_', run_label, '.fig'])
% % Read scalar field from VTK (handles line-wrapped ASCII data automatically)
% fid = fopen('atp_cyto_final.vtk', 'r');
% n_points = 0;
% tline = fgetl(fid);
% while ischar(tline)
%     if startsWith(tline, 'POINT_DATA')
%         parts = strsplit(tline);
%         n_points = str2double(parts{2});
%     end
%     if startsWith(tline, 'LOOKUP_TABLE')
%         break;
%     end
%     tline = fgetl(fid);
% end
% atp_vals = fscanf(fid, '%f', n_points);
% fclose(fid);
% 
% % Producer mask must be in the same flattened CRU order as the VTK point data
% mask = readmatrix('producer_mask.txt', 'FileType', 'text', 'Delimiter', '\t');
% mask = mask(:);
% 
% prod_cols    = find(mask == 1);
% nonprod_cols = find(mask == 0);
% 
% atp_prod    = atp_vals(prod_cols);
% atp_nonprod = atp_vals(nonprod_cols);
% 
% atp_prod_mean    = mean(atp_prod);
% atp_prod_std     = std(atp_prod);
% atp_nonprod_mean = mean(atp_nonprod);
% atp_nonprod_std  = std(atp_nonprod);
% 
% fprintf('Producer ATP:     mean = %.4f uM, std = %.4f uM (n=%d)\n', atp_prod_mean, atp_prod_std, numel(atp_prod));
% fprintf('Non-producer ATP: mean = %.4f uM, std = %.4f uM (n=%d)\n', atp_nonprod_mean, atp_nonprod_std, numel(atp_nonprod));
% fprintf('Difference (Producer - Non-producer): %.4f uM\n', atp_prod_mean - atp_nonprod_mean);
% 
% means = [atp_prod_mean, atp_nonprod_mean];
% stds  = [atp_prod_std, atp_nonprod_std];
% 
% figure('Name', 'Producer vs Non-Producer ATP', 'Position', [1050 1300 500 400], 'Color', 'w');
% b = bar(means, 'FaceColor', 'flat');
% b.CData(1,:) = [0.00 0.45 0.74];
% b.CData(2,:) = [0.85 0.33 0.10];
% hold on;
% errorbar(1:2, means, stds, 'k', 'LineStyle', 'none', 'LineWidth', 1.2);
% hold off;
% 
% set(gca, 'XTick', 1:2, 'XTickLabel', {'Producer', 'Non-producer'});
% ylabel('ATP (\muM)');
% title('Producer vs. Non-Producer ATP (Final Snapshot)');
% %savefig(['ATP_Bar_Final_', run_label, '.fig']);
% 
% atp_vals = read_vtk_scalars('atp_cyto_final.vtk');
% status   = read_vtk_scalars('producer_status_map.vtk');
% groups = {'Non-producer','Producer (unimpaired)','Producer (impaired)'};
% means = zeros(1,3); stds = zeros(1,3);
% for g = 0:2
%     vals = atp_vals(status == g);
%     means(g+1) = mean(vals);
%     stds(g+1)  = std(vals);
% end
% 
% figure;
% bar(means); hold on;
% errorbar(1:3, means, stds, 'k.', 'LineWidth', 1);
% set(gca, 'XTickLabel', groups);
% ylabel('ATP (\muM)');
% title('Final ATP by CRU category');
% %savefig(['atp_bar_by_category', run_label, '.fig']);



%% ATP-dependent current traces (ikatp_inak_serca_trace.txt)
% Columns: time  ATP_ave  p_kATP  ikATP  inak  Iup_avg
trace_data = readtable('ikatp_inak_serca_trace.txt', 'FileType', 'text', 'Delimiter', '\t');
t_trace  = trace_data.time;
p_katp   = trace_data.p_kATP;
ikatp    = trace_data.ikATP;
inak_tr  = trace_data.inak;
iup_avg  = trace_data.Iup_avg;

figure('Name', 'ATP-Dependent Currents', 'Position', [100 100 900 900], 'Color', 'w');
subplot(3,1,1);
plot(t_trace, ikatp, 'LineWidth', 1.2, 'Color', [0.85 0.33 0.10]);
ylabel('I_{KATP} (pA)');
title('ATP-Dependent Current Traces');

subplot(3,1,2);
plot(t_trace, inak_tr, 'LineWidth', 1.2, 'Color', [0.00 0.45 0.74]);
ylabel('I_{NaK} (pA)');

subplot(3,1,3);
plot(t_trace, iup_avg, 'LineWidth', 1.2, 'Color', [0.47 0.67 0.19]);
ylabel('I_{up,avg} (\muM/ms)');
xlabel('Time (ms)');
savefig(['ATPdep_currents_', run_label, '.fig'])
 
% %Line graph by prod status 
% data = readmatrix('atp_group_trace.txt', 'FileType', 'text', 'Delimiter', '\t', 'NumHeaderLines', 1);
% t = data(:,1);
% figure;
% plot(t, data(:,2), 'DisplayName', 'Non-producer'); hold on;
% plot(t, data(:,3), 'DisplayName', 'Producer (unimpaired)');
% plot(t, data(:,4), 'DisplayName', 'Producer (impaired)');
% xlabel('Time (ms)'); ylabel('ATP (\muM)');
% legend show;
% title('ATP by CRU category over time');
% %savefig(['atp_line_by_category', run_label, '.fig']);
% 
% %Bar graph with 3 prod. categories 
% atp_vals = read_vtk_scalars('atp_cyto_final.vtk');
% status   = read_vtk_scalars('producer_status_map.vtk');
% 
% groups = {'Non-producer','Producer (unimpaired)','Producer (impaired)'};
% means = zeros(1,3); stds = zeros(1,3);
% for g = 0:2
%     vals = atp_vals(status == g);
%     means(g+1) = mean(vals);
%     stds(g+1)  = std(vals);
% end
% 
% figure;
% bar(means); hold on;
% errorbar(1:3, means, stds, 'k.', 'LineWidth', 1);
% set(gca, 'XTickLabel', groups);
% ylabel('ATP (\muM)');
% title('Final ATP by CRU category');
% savefig(['atp_bar_by_3category', run_label, '.fig']);
% 
% %Strip Plot of Final ATP
% status   = read_vtk_scalars('producer_status_map.vtk');
% groups = categorical(status, [0 1 2], {'Non-producer','Producer (unimpaired)','Producer (impaired)'});
% figure;
% swarmchart(groups, atp_vals, 15, 'filled');
% ylabel('ATP (\muM)');
% title('ATP distribution by CRU category (final snapshot)');
% hold on;% same 20 s window
% cats = categories(groups);
% for i = 1:numel(cats)
%     m = mean(atp_vals(groups == cats{i}));
%     plot([i-0.25 i+0.25], [m m], 'k-', 'LineWidth', 2);
% end
% %savefig(['atp_strip_by_3category', run_label, '.fig']);
%% ===== Settings used by every block below =====
win_ms       = 20000;   % analysis window: last 20 s (use 30000 for slowed runs)
mito_per_CRU = 1/4;     % 4:1 run -> 1/4 ; 8:1 run -> 1/8

w  = t       >= t(end)       - win_ms;     % mito_summary window
wt = t_trace >= t_trace(end) - win_ms;     % current-trace window
wv = t_ci    >= t_ci(end)    - win_ms;     % ci.txt window

row = struct();
row.label        = {run_label};
row.win_s        = win_ms/1000;
row.mito_per_CRU = mito_per_CRU;

%% Beat-based cytosolic biomarkers (NaN if the cell is not beating)
cat_names = {'diastolic_Ca','systolic_Ca','CaT_amplitude','CaD50'};
ap_names  = {'MDP','AP_amplitude','APD90'};     % from results.ap_markers
beat_ok = true;
try
    results = biomarkers('ci.txt', 'ikatp_inak_serca_trace.txt', run_label, win_ms, false);
    if isempty(results.BCL_vm) || all(isnan(results.BCL_vm)), beat_ok = false; end
catch ME
    fprintf('biomarkers() failed (%s): treating run as non-beating\n', ME.message);
    beat_ok = false;
end

if beat_ok
    row.BCL_vm_mean   = mean(results.BCL_vm, 'omitnan');
    row.BCL_vm_sd     = std(results.BCL_vm, 'omitnan');
    row.BCL_vm_median = median(results.BCL_vm, 'omitnan');
    row.CaT_CL_mean   = mean(results.CaT_CL, 'omitnan');
    row.CaT_CL_sd     = std(results.CaT_CL, 'omitnan');
    for i = 1:numel(cat_names)
        vals = results.cat_markers.(cat_names{i});
        row.([cat_names{i} '_mean']) = mean(vals, 'omitnan');
        row.([cat_names{i} '_sd'])   = std(vals, 'omitnan');
    end
    for i = 1:numel(ap_names)
        if isfield(results, 'ap_markers') && isfield(results.ap_markers, ap_names{i})
            vals = results.ap_markers.(ap_names{i});
            row.([ap_names{i} '_mean']) = mean(vals, 'omitnan');
            row.([ap_names{i} '_sd'])   = std(vals, 'omitnan');
        else
            row.([ap_names{i} '_mean']) = NaN;   % older biomarkers.m without AP markers
            row.([ap_names{i} '_sd'])   = NaN;
        end
    end
else
    [row.BCL_vm_mean, row.BCL_vm_sd, row.BCL_vm_median, ...
     row.CaT_CL_mean, row.CaT_CL_sd] = deal(NaN);
    for i = 1:numel(cat_names)
        row.([cat_names{i} '_mean']) = NaN;
        row.([cat_names{i} '_sd'])   = NaN;
    end
    for i = 1:numel(ap_names)
        row.([ap_names{i} '_mean']) = NaN;
        row.([ap_names{i} '_sd'])   = NaN;
    end
end

%% ATP-dependent currents: window means (valid whether beating or not)
cur_names = {'ATP_ave','p_kATP','ikATP','inak','Iup_avg'};
for i = 1:numel(cur_names)
    vals = trace_data.(cur_names{i})(wt);
    row.([cur_names{i} '_mean']) = mean(vals, 'omitnan');
    row.([cur_names{i} '_sd'])   = std(vals, 'omitnan');
end
row.ATP_ave_max = max(trace_data.ATP_ave(wt), [], 'omitnan');
row.ATP_ave_min = min(trace_data.ATP_ave(wt), [], 'omitnan');

%% Cycle length from upstrokes (threshold crossing and max dV/dt)
try
    cl = cl_upstroke(t_ci, vm_ci, win_ms);
catch
    cl = struct('n_beats', 0, 'CL_thr', [], 'CL_dvdt', [], 't_beats', []);
end
row.n_beats        = cl.n_beats;
row.CL_thr_mean    = mean(cl.CL_thr,    'omitnan');
row.CL_thr_median  = median(cl.CL_thr,  'omitnan');
row.CL_thr_sd      = std(cl.CL_thr,     'omitnan');
row.CL_dvdt_mean   = mean(cl.CL_dvdt,   'omitnan');
row.CL_dvdt_median = median(cl.CL_dvdt, 'omitnan');
row.CL_dvdt_sd     = std(cl.CL_dvdt,    'omitnan');
rate_Hz = 1000 ./ cl.CL_thr ;
row.rate_Hz_mean  = mean(rate_Hz, 'omitnan');
row.rate_Hz_sd    = std(rate_Hz, 'omitnan');
row.rate_bpm_mean = 60 * row.rate_Hz_mean;
row.rate_bpm_sd   = 60 * row.rate_Hz_sd;

%% Vm window statistics and cell state
Vm_w = vm_ci(wv);
row.Vm_mean  = mean(Vm_w);
row.Vm_sd    = std(Vm_w);
row.Vm_min   = min(Vm_w);
row.Vm_max   = max(Vm_w);
Vm_range     = row.Vm_max - row.Vm_min;

if cl.n_beats >= 3 && Vm_range > 30
    if row.CL_thr_sd / row.CL_thr_mean < 0.10
        state = "beating";
    else
        state = "irregular";
    end
elseif Vm_range < 5
    state = "arrested";
else
    state = "subthreshold";             % oscillating, no full APs
end
row.state = {char(state)};
fprintf('\nState: %s | Vm mean %.2f mV (min %.2f, max %.2f)\n', ...
        state, row.Vm_mean, row.Vm_min, row.Vm_max);

% No real beats -> beat-based values are meaningless, blank them
if state == "arrested" || state == "subthreshold"
    blank = [{'BCL_vm_mean','BCL_vm_sd','BCL_vm_median','CaT_CL_mean','CaT_CL_sd', ...
              'CL_thr_mean','CL_thr_median','CL_thr_sd', ...
              'CL_dvdt_mean','CL_dvdt_median','CL_dvdt_sd', 'rate_Hz_mean','rate_Hz_sd','rate_bpm_mean','rate_bpm_sd'}, ...
             strcat(cat_names, '_mean'), strcat(cat_names, '_sd'), ...
             strcat(ap_names,  '_mean'), strcat(ap_names,  '_sd')];
    for i = 1:numel(blank), row.(blank{i}) = NaN; end
end

%% Window means from mito_summary (whole-cell averages; valid in any state)
row.avg_cai_mean     = mean(avg_cai(w));
row.avg_cai_sd       = std(avg_cai(w));
row.avg_cp_prod_mean = mean(avg_cp(w));
row.avg_ca_mito_mean = mean(ca_mito(w));
row.avg_ca_mito_sd   = std(ca_mito(w));
row.avg_psi_mean     = mean(psi_mito(w));
row.avg_psi_sd       = std(psi_mito(w));
row.avg_ATP_mean     = mean(atp(w));
row.avg_ADP_mean     = mean(adp(w));
p = polyfit(t(w)/1000, atp(w), 1);          % ATP drift over the window
row.ATP_slope_uM_per_s = p(1);
row.J_uni_mean       = mean(data.avg_J_uni(w));
row.jNCX_m_mean      = mean(data.avg_jNCX_m(w));

%% Mito Ca shuttle vs SERCA
vp     = 0.00126 * 2;        % cleft volume per CRU, um^3
vi     = 1.5 * 0.5;          % cytosolic volume per CRU, um^3
BCL    = row.BCL_vm_mean;    % NaN when not beating

Jflux_cp   = mean(data.avg_Jflux_cp(w));
NCXflux_ci = mean(data.avg_NCXflux_ci(w));
Iup        = mean(iup_avg(wt));

% Rates (zmol per second) - defined in any state
mcu_rate   = Jflux_cp   * vp * 1000;
nclx_rate  = NCXflux_ci * vi * 1000;
serca_rate = Iup        * vi * 1000;

% Ratios do not depend on BCL, so they survive arrest
shuttle_frac = (mcu_rate * mito_per_CRU) / serca_rate;

% Per-beat amounts (NaN when not beating)
mcu_per_mito  = mcu_rate   * BCL / 1000;
nclx_per_mito = nclx_rate  * BCL / 1000;
serca_per_CRU = serca_rate * BCL / 1000;

fprintf('\n--- Mito Ca shuttle (last %.0f s) ---\n', win_ms/1000);
fprintf('MCU uptake per mito  : %.3f zmol/s | %.3f zmol/beat\n', mcu_rate,   mcu_per_mito);
fprintf('NCLX return per mito : %.3f zmol/s | %.3f zmol/beat\n', nclx_rate,  nclx_per_mito);
fprintf('SERCA uptake per CRU : %.3f zmol/s | %.3f zmol/beat\n', serca_rate, serca_per_CRU);
fprintf('Mito shuttle / SERCA (whole cell): %.2f %%\n', 100*shuttle_frac);

row.Jflux_cp_mean       = Jflux_cp;
row.NCXflux_ci_mean     = NCXflux_ci;
row.mcu_rate_zmol_s     = mcu_rate;
row.nclx_rate_zmol_s    = nclx_rate;
row.serca_rate_zmol_s   = serca_rate;
row.mcu_per_mito_zmol   = mcu_per_mito;
row.nclx_per_mito_zmol  = nclx_per_mito;
row.uptake_return_ratio = mcu_rate / nclx_rate;   % ~1 at steady state
row.serca_per_CRU_zmol  = serca_per_CRU;
row.shuttle_frac_pct    = 100*shuttle_frac;

%% Mito biomarkers per beat (NaN when there are no beats)
tb = cl.t_beats;  nb = max(numel(tb) - 1, 0);
caDia = nan(nb,1); caPeak = nan(nb,1); psiRest = nan(nb,1); psiNadir = nan(nb,1);
if state == "beating" || state == "irregular"
    for k = 1:nb
        m = t >= tb(k) & t < tb(k+1);
        if nnz(m) < 5, continue; end
        caDia(k)    = min(ca_mito(m));     caPeak(k)   = max(ca_mito(m));
        psiRest(k)  = max(psi_mito(m));    psiNadir(k) = min(psi_mito(m));
    end
end
caAmp = caPeak - caDia;   dPsi = psiRest - psiNadir;

s = @(x) [mean(x,'omitnan') std(x,'omitnan')];
fprintf('  Mito Ca (window mean): %.1f ± %.1f nM\n', 1000*mean(ca_mito(w)), 1000*std(ca_mito(w)));
fprintf('  Psi (window mean)    : %.2f ± %.2f mV\n', mean(psi_mito(w)), std(psi_mito(w)));
fprintf('  Mito Ca diastolic : %.1f ± %.1f nM\n', 1000*s(caDia));
fprintf('  Mito Ca peak      : %.1f ± %.1f nM\n', 1000*s(caPeak));
fprintf('  Mito CaT amplitude: %.1f ± %.1f nM\n', 1000*s(caAmp));
fprintf('  Psi resting (max) : %.2f ± %.2f mV\n', s(psiRest));
fprintf('  Psi nadir (min)   : %.2f ± %.2f mV\n', s(psiNadir));
fprintf('  dPsi per beat     : %.2f ± %.2f mV  (n=%d beats)\n', s(dPsi), nnz(~isnan(dPsi)));

%% Mito by CRU type (0 interior, 1 near-boundary, 2 boundary)
if isfile('mito_by_type.txt')
    ty = readtable('mito_by_type.txt', 'FileType', 'text', 'Delimiter', '\t');
    tt = ty.time;  names = {'Interior','Near-boundary','Boundary'};
    cols = [0.00 0.45 0.74; 0.93 0.69 0.13; 0.85 0.33 0.10];

    figure('Name','Mito by CRU type','Position',[1050 100 900 900],'Color','w');
    vars = {'cp','ca_mito','psi'};  ylabs = {'Cleft Ca (\muM)','Mito Ca (\muM)','\Psi_{mito} (mV)'};
    for p = 1:3
        ax(p) = subplot(3,1,p); hold on;
        for k = 0:2, plot(tt, ty.(sprintf('%s_%d', vars{p}, k)), 'Color', cols(k+1,:), 'LineWidth', 1.2); end
        ylabel(ylabs{p});
    end
    legend(ax(1), names, 'Location', 'best');  xlabel('Time (ms)');  linkaxes(ax, 'x');
    savefig(['mito_by_type_', run_label, '.fig'])

    fprintf('\n--- Per CRU type (last %.0f s) ---\n', win_ms/1000);
    for k = 0:2
        cm = ty.(sprintf('ca_mito_%d', k));  ps = ty.(sprintf('psi_%d', k));  cp = ty.(sprintf('cp_%d', k));
        amp = nan(nb,1); dps = nan(nb,1); cpk = nan(nb,1);
        if state == "beating" || state == "irregular"
            for b = 1:nb
                m = tt >= tb(b) & tt < tb(b+1);
                if nnz(m) < 5, continue; end
                amp(b) = max(cm(m)) - min(cm(m));  dps(b) = max(ps(m)) - min(ps(m));  cpk(b) = max(cp(m));
            end
        end
        fprintf('%-14s mito swing %6.2f ± %.2f nM | dPsi %.3f ± %.3f mV | peak cleft %6.2f ± %.2f uM\n', ...
            names{k+1}, 1000*s(amp), s(dps), s(cpk));
        tag = sprintf('type%d', k);
        row.(['mitoCaT_amp_nM_' tag]) = 1000*mean(amp,'omitnan');  row.(['mitoCaT_amp_sd_' tag]) = 1000*std(amp,'omitnan');
        row.(['dPsi_mV_' tag])        = mean(dps,'omitnan');       row.(['dPsi_sd_' tag])        = std(dps,'omitnan');
        row.(['cleft_peak_uM_' tag])  = mean(cpk,'omitnan');       row.(['cleft_peak_sd_' tag])  = std(cpk,'omitnan');
    end
end

row.mitoCa_dia_nM  = 1000*mean(caDia,'omitnan');   row.mitoCa_dia_sd  = 1000*std(caDia,'omitnan');
row.mitoCa_peak_nM = 1000*mean(caPeak,'omitnan');  row.mitoCa_peak_sd = 1000*std(caPeak,'omitnan');
row.mitoCaT_amp_nM = 1000*mean(caAmp,'omitnan');   row.mitoCaT_amp_sd = 1000*std(caAmp,'omitnan');
row.psi_rest_mV    = mean(psiRest,'omitnan');      row.psi_rest_sd    = std(psiRest,'omitnan');
row.psi_nadir_mV   = mean(psiNadir,'omitnan');     row.psi_nadir_sd   = std(psiNadir,'omitnan');
row.dPsi_mV        = mean(dPsi,'omitnan');         row.dPsi_sd        = std(dPsi,'omitnan');

%% Write everything (last step)
T = struct2table(row);
writetable(T, sprintf('biomarkers_%s.csv', run_label));
fprintf('Saved biomarker summary: biomarkers_%s.csv (%d columns, state = %s)\n', ...
        run_label, width(T), state);

function vals = read_vtk_scalars(filename)
fid = fopen(filename, 'r');
tline = fgetl(fid);
n_points = 0;
while ischar(tline)
    if startsWith(tline, 'POINT_DATA')
        parts = strsplit(tline);
        n_points = str2double(parts{2});
    end
    if startsWith(tline, 'LOOKUP_TABLE'), break; end
    tline = fgetl(fid);
end
vals = fscanf(fid, '%f', n_points);
fclose(fid);
end