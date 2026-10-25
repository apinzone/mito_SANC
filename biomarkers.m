function results = biomarkers(ci_path, currents_path, label, window_ms, save_mat)
%ANALYZE_MODEL  Compute BCL, CaT amplitude, diastolic/systolic Ca, CaD50 (and,
%   if a currents trace is supplied, current summaries) for ONE model run,
%   restricted to the last WINDOW_MS of the simulation.
%
%   results = analyze_model(ci_path, currents_path, label, window_ms, save_mat)
%
%   ci_path       - path to this model's ci.txt (time, Vm, Cai). Use 'ci.txt'
%                   if running from inside the model's own output directory.
%   currents_path - path to this model's currents trace:
%                     'ikatp_inak_serca_trace.txt' for the updated model
%                     'inak_iup_trace.txt'         for the backup model
%                   pass '' to skip.
%   label         - short string identifying this run, e.g. 'Backup' or
%                   'Updated' (used in plot legends / printed output / the
%                   saved .mat filename).
%   window_ms     - analysis window: last N ms of the sim (default 100000)
%   save_mat      - true/false, save results to '<label>_results.mat' in the
%                   current directory (default true) so both models' results
%                   can be loaded together later. Default true.
%
%   Run this from inside each model's own output folder:
%     analyze_model('ci.txt', 'inak_iup_trace.txt', 'Backup', 100000);
%     analyze_model('ci.txt', 'ikatp_inak_serca_trace.txt', 'Updated', 100000);
%
%   Then, from wherever the two .mat files end up:
%     compare_models('Backup_results.mat', 'Updated_results.mat')

if nargin < 2, currents_path = ''; end
if nargin < 3 || isempty(label), label = 'Model'; end
if nargin < 4 || isempty(window_ms), window_ms = 100000; end
if nargin < 5 || isempty(save_mat), save_mat = true; end

is_updated = ~isempty(currents_path) && contains(currents_path, 'ikatp');

fprintf('\n=== %s (last %.0f s) ===\n', label, window_ms/1000);

data = readmatrix(ci_path, 'FileType', 'text', 'Delimiter', ' ');
t = data(:,1); vm = data(:,2); cai = data(:,3);

mask = t >= (t(end) - window_ms);
t = t(mask); vm = vm(mask); cai = cai(mask);

% --- BCL from Vm peaks ---
[peaks_vm, peak_t_vm] = find_beat_peaks(t, vm);
cl_vm = diff(peak_t_vm);
fprintf('  Beats detected (Vm): %d\n', numel(peak_t_vm));
summarize(cl_vm, 'BCL (Vm peak-to-peak, ms)');

% --- Ca transient peaks + per-beat biomarkers ---
[peaks_ca, peak_t_ca] = find_beat_peaks(t, cai);
cl_ca = diff(peak_t_ca);
summarize(cl_ca, 'CaT CL (Cai peak-to-peak, ms)');

cat_markers = compute_cat_biomarkers(t, cai, peaks_ca);
fn = fieldnames(cat_markers);
for i = 1:numel(fn)
    summarize(cat_markers.(fn{i}), fn{i});
end

% --- AP biomarkers (MDP, AP amplitude, APD90) ---
ap_markers = compute_ap_biomarkers(t, vm, peaks_vm);
fn = fieldnames(ap_markers);
for i = 1:numel(fn)
    summarize(ap_markers.(fn{i}), fn{i});
end

results.t = t; results.vm = vm; results.cai = cai;
results.BCL_vm = cl_vm; results.CaT_CL = cl_ca;
results.cat_markers = cat_markers;
results.ap_markers  = ap_markers;
results.label = label;
results.currents = struct();

% --- Currents (optional) ---
if ~isempty(currents_path) && isfile(currents_path)
    cdata = readtable(currents_path, 'FileType', 'text', 'Delimiter', '\t');
    ct = cdata.time;
    cmask = ct >= (ct(end) - window_ms);
    cdata = cdata(cmask, :);
    results.currents.time = cdata.time;   % masked time vector, needed for plotting
    cols = {'ATP_ave','p_kATP','ikATP','inak','Iup_avg'};
    for i = 1:numel(cols)
        c = cols{i};
        if ismember(c, cdata.Properties.VariableNames)
            vals = cdata.(c);
            summarize(vals, c);
            results.currents.(c) = vals;
        end
    end
end

if save_mat
    outname = sprintf('%s_results.mat', matlab.lang.makeValidName(label));
    save(outname, 'results');
    fprintf('\nSaved results: %s\n', outname);
end

end

%% ------------------------------------------------------------------------
function [peaks_idx, peak_times] = find_beat_peaks(t, signal, min_cl_guess_ms, prominence_frac, height_frac)
if nargin < 3 || isempty(min_cl_guess_ms), min_cl_guess_ms = 100; end
if nargin < 4 || isempty(prominence_frac), prominence_frac = 0.3; end
if nargin < 5 || isempty(height_frac),     height_frac = 0.5;     end

sig_min = min(signal); sig_max = max(signal);
height_val     = sig_min + height_frac * (sig_max - sig_min);
prominence_val = prominence_frac * (sig_max - sig_min);

dt = median(diff(t));
distance_samples = max(1, round(min_cl_guess_ms / dt));

[~, peaks_idx] = findpeaks(signal, 'MinPeakHeight', height_val, ...
    'MinPeakProminence', prominence_val, 'MinPeakDistance', distance_samples);

peak_times = t(peaks_idx);
end

%% ------------------------------------------------------------------------
function markers = compute_cat_biomarkers(t, cai, peaks_idx, decay_pct)
%COMPUTE_CAT_BIOMARKERS  Per-beat diastolic/systolic Ca, amplitude, CaD(decay_pct).
%   CaD50 (decay_pct = 0.5) is the Ca TRANSIENT DURATION at 50% amplitude: the
%   time from the up-crossing of (diastolic + 0.5*amplitude) on the rising phase
%   to the down-crossing of the same level on the falling phase (full width at
%   half maximum). Both crossings are linearly interpolated.
%   NaN for any beat whose crossings are incomplete (record starts/ends mid-beat).
%   (Changed 2026-10-23 from the previous peak -> 50% decay definition.)
if nargin < 4 || isempty(decay_pct), decay_pct = 0.5; end

n = numel(peaks_idx);
diastolic = nan(max(n-1,0), 1);
systolic  = nan(max(n-1,0), 1);
amplitude = nan(max(n-1,0), 1);
cad       = nan(max(n-1,0), 1);

for i = 1:(n-1)                      % only beats with a following peak -> complete beats
    p0 = peaks_idx(i); p1 = peaks_idx(i+1);
    seg = cai(p0:p1);

    diastolic(i) = min(seg);
    systolic(i)  = cai(p0);
    amplitude(i) = systolic(i) - diastolic(i);

    cad(i) = beat_width(t, cai, peaks_idx, i, decay_pct);
end

label = sprintf('CaD%d', round(decay_pct*100));
markers.diastolic_Ca  = diastolic;
markers.systolic_Ca   = systolic;
markers.CaT_amplitude = amplitude;
markers.(label)       = cad;
end

%% ------------------------------------------------------------------------
function markers = compute_ap_biomarkers(t, vm, peaks_idx, repol_pct)
%COMPUTE_AP_BIOMARKERS  Per-beat MDP, AP peak/amplitude and APD(repol_pct).
%   APD90 (repol_pct = 0.9) runs from the upstroke (time of maximum dV/dt) to the
%   down-crossing of MDP + (1-0.9)*(peak - MDP) on the repolarizing phase, with
%   the crossing linearly interpolated.
%   NOTE: a symmetric up/down-crossing width is NOT used here. In a spontaneously
%   firing cell, diastolic depolarization carries Vm above the 90% repolarization
%   level long before the upstroke, so a width measure would include most of
%   diastole. MDP is the trough that follows each AP.
%   NaN for any beat whose crossings are incomplete (record ends mid-beat).
if nargin < 4 || isempty(repol_pct), repol_pct = 0.9; end

n = numel(peaks_idx);
mdp   = nan(max(n-1,0), 1);
vpeak = nan(max(n-1,0), 1);
apa   = nan(max(n-1,0), 1);
apd   = nan(max(n-1,0), 1);

for i = 1:(n-1)
    p0 = peaks_idx(i); p1 = peaks_idx(i+1);
    seg = vm(p0:p1); seg_t = t(p0:p1);

    mdp(i)   = min(seg);
    vpeak(i) = vm(p0);
    apa(i)   = vpeak(i) - mdp(i);

    % upstroke: max dV/dt between the preceding trough and this peak
    if i > 1, q0 = peaks_idx(i-1); else, q0 = max(1, p0 - (p1 - p0)); end
    pre = vm(q0:p0); pre_t = t(q0:p0);
    [~, kmin] = min(pre);
    rise = pre(kmin:end); rise_t = pre_t(kmin:end);
    if numel(rise) < 3, continue; end
    dv = diff(rise) ./ diff(rise_t);
    [~, ku] = max(dv);
    t_up = 0.5*(rise_t(ku) + rise_t(ku+1));

    % 90% repolarization level, down-crossing on the falling phase
    lvl = mdp(i) + (1 - repol_pct) * apa(i);
    kd  = find(seg(1:end-1) > lvl & seg(2:end) <= lvl, 1, 'first');
    if isempty(kd), continue; end                       % never repolarizes -> NaN
    t_dn = seg_t(kd) + (lvl - seg(kd)) / (seg(kd+1) - seg(kd)) * (seg_t(kd+1) - seg_t(kd));

    apd(i) = t_dn - t_up;
end

label = sprintf('APD%d', round(repol_pct*100));
markers.MDP          = mdp;
markers.AP_peak      = vpeak;
markers.AP_amplitude = apa;
markers.(label)      = apd;
end

%% ------------------------------------------------------------------------
function w = beat_width(t, y, peaks_idx, i, pct)
%BEAT_WIDTH  Duration of beat i at the level base + (1-pct)*(peak - base).
%   base is the minimum between the previous peak (or one cycle of lead-in when
%   beat i is the first) and this peak, so the rising phase is always inside the
%   search window. Both crossings are linearly interpolated; NaN if either is
%   missing, which also covers a record that starts or ends mid-beat.
p0 = peaks_idx(i); p1 = peaks_idx(i+1);
if i > 1, q0 = peaks_idx(i-1); else, q0 = max(1, p0 - (p1 - p0)); end

pre = y(q0:p0); pre_t = t(q0:p0);
[base, kmin] = min(pre);
lvl = base + (1 - pct) * (y(p0) - base);

rise = pre(kmin:end);  rise_t = pre_t(kmin:end);
fall = y(p0:p1);       fall_t = t(p0:p1);

ku = find(rise(1:end-1) <  lvl & rise(2:end) >= lvl, 1, 'last');    % upstroke crossing
kd = find(fall(1:end-1) >  lvl & fall(2:end) <= lvl, 1, 'first');   % downstroke crossing

w = NaN;
if ~isempty(ku) && ~isempty(kd)
    t_up = rise_t(ku) + (lvl - rise(ku)) / (rise(ku+1) - rise(ku)) * (rise_t(ku+1) - rise_t(ku));
    t_dn = fall_t(kd) + (lvl - fall(kd)) / (fall(kd+1) - fall(kd)) * (fall_t(kd+1) - fall_t(kd));
    w = t_dn - t_up;
end
end

%% ------------------------------------------------------------------------
function stats = summarize(arr, name)
arr = arr(~isnan(arr));
if isempty(arr)
    fprintf('  %s: no valid beats\n', name);
    stats = [];
    return;
end
stats.mean = mean(arr); stats.sd = std(arr); stats.median = median(arr); stats.n = numel(arr);
fprintf('  %s: mean=%.4f  SD=%.4f  median=%.4f  n=%d\n', ...
    name, stats.mean, stats.sd, stats.median, stats.n);
end