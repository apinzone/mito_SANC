function out = cl_upstroke(t, v, win_ms, v_thr)
%CL_UPSTROKE Cycle length from AP upstrokes over the last win_ms of a run.
%   out.CL_dvdt : intervals between times of max dV/dt (one per beat)
%   out.CL_thr  : intervals between upward crossings of v_thr (interpolated)
%   t in ms, v in mV. Default v_thr = -20 mV.
if nargin < 4, v_thr = -20; end
t = t(:); v = v(:);
w = t >= t(end) - win_ms;
t = t(w); v = v(w);

% Upward threshold crossings, with linear interpolation for sub-sample timing
up = find(v(1:end-1) < v_thr & v(2:end) >= v_thr);
t_thr = t(up) + (v_thr - v(up)) .* (t(up+1) - t(up)) ./ (v(up+1) - v(up));

% Max dV/dt within a short window around each crossing (the upstroke)
dvdt = gradient(v, t);
t_dvdt = nan(size(up));
for k = 1:numel(up)
    i0 = max(1, up(k) - 5);  i1 = min(numel(v), up(k) + 10);
    [~, im] = max(dvdt(i0:i1));
    t_dvdt(k) = t(i0 + im - 1);
end

out.CL_thr  = diff(t_thr);
out.CL_dvdt = diff(t_dvdt);
out.n_beats = numel(up);
f = @(x) [mean(x,'omitnan') median(x,'omitnan') std(x,'omitnan')];
s1 = f(out.CL_thr);  s2 = f(out.CL_dvdt);
fprintf('  CL (%g mV crossing): mean=%.2f  median=%.2f  SD=%.2f ms  n=%d\n', v_thr, s1, numel(out.CL_thr));
fprintf('  CL (max dV/dt)     : mean=%.2f  median=%.2f  SD=%.2f ms  n=%d\n', s2, numel(out.CL_dvdt));
end