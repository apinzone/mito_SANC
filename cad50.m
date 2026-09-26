function cad50 = compute_cad50(ci_path, window_ms, decay_pct)
%COMPUTE_CAD50  Quick standalone CaD50 calculation from a single ci.txt.
%   Same corrected logic as analyze_model.m: full width of the Ca transient
%   at the decay_pct level (rising crossing before the peak to falling
%   crossing after the peak), not just peak-to-decay.
%
%   cad50 = compute_cad50('ci.txt')
%   cad50 = compute_cad50('ci.txt', 100000)      % last 100s only
%   cad50 = compute_cad50('ci.txt', 100000, 0.9) % CaD90 instead of CaD50

if nargin < 2 || isempty(window_ms), window_ms = 100000; end
if nargin < 3 || isempty(decay_pct), decay_pct = 0.5; end

data = readmatrix(ci_path, 'FileType', 'text', 'Delimiter', ' ');
t = data(:,1); cai = data(:,3);

mask = t >= (t(end) - window_ms);
t = t(mask); cai = cai(mask);

% --- peak detection (same thresholds as analyze_model.m) ---
sig_min = min(cai); sig_max = max(cai);
height_val = sig_min + 0.5 * (sig_max - sig_min);
prominence_val = 0.3 * (sig_max - sig_min);
dt = median(diff(t));
distance_samples = max(1, round(100 / dt));

[~, peaks_idx] = findpeaks(cai, 'MinPeakHeight', height_val, ...
    'MinPeakProminence', prominence_val, 'MinPeakDistance', distance_samples);

n = numel(peaks_idx);
cad = nan(n-2, 1);

for k = 1:(n-2)
    i = k + 1;
    p_prev = peaks_idx(i-1); p0 = peaks_idx(i); p_next = peaks_idx(i+1);

    pre_seg = cai(p_prev:p0);
    pre_t   = t(p_prev:p0);
    [diastolic, trough_local] = min(pre_seg);
    systolic  = cai(p0);
    amplitude = systolic - diastolic;
    target = diastolic + (1 - decay_pct) * amplitude;

    rise_seg = pre_seg(trough_local:end);
    rise_t   = pre_t(trough_local:end);
    above = find(rise_seg > target, 1, 'first');
    if isempty(above), continue; end
    t_rise = rise_t(above);

    fall_seg = cai(p0:p_next);
    fall_t   = t(p0:p_next);
    below = find(fall_seg <= target, 1, 'first');
    if isempty(below), continue; end
    t_fall = fall_t(below);

    cad(k) = t_fall - t_rise;
end

cad50 = cad(~isnan(cad));
label = sprintf('CaD%d', round(decay_pct*100));
fprintf('%s: mean=%.4f  SD=%.4f  median=%.4f  n=%d (of %d beats)\n', ...
    label, mean(cad50), std(cad50), median(cad50), numel(cad50), n-2);

end