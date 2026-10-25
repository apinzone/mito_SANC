%% make_AP_templates.m
% Builds AP-clamp waveforms from a spontaneous control run and overlays
% them on the control beats.
%
% Input : vm_hr.txt  (t [ms], Vm [mV]) from an UNBLOCKED spontaneous run,
%         last ~2 s written every 10 steps (0.1 ms)
% Output: pool_AP/AP_<CL>.txt, one Vm value per line, sampled every dt_out,
%         exactly round(CL/dt_out) lines, starting at the MDP
%
% Template = [diastole (MDP -> split)] + [AP (split -> next MDP), shape fixed]
% split = fixed time before the -20 mV upstroke (upstroke-locked average)
%
% Two ways to fill the diastole up to the target CL:
%   hold_mode = false : stretch the native diastolic ramp uniformly
%                       (use for physiological CLs, ~106-350 ms)
%   hold_mode = true  : hold Vm at the MDP, then play the NATIVE diastolic
%                       ramp unchanged before the AP
%                       (use for very long CLs, e.g. 2000-5000 ms, where a
%                        stretched ramp would be unrealistically slow)

clear; clc;

%% Settings
infile    = 'vm_hr_CONTROL.txt';
CL_list   = [10000];         % target cycle lengths (ms)
hold_mode = false;           % see header
dt_out    = 0.01;           % must equal dt in SAN_pace.cpp
v_thr     = -20;            % upstroke detection threshold (mV)
dvdt_TO   = 0.5;            % takeoff criterion (mV/ms)
outdir    = 'pool_AP';

%% Load
d  = readmatrix(infile);
t  = d(:,1);  V = d(:,2);
dt_in = median(diff(t));
dVdt  = gradient(V, dt_in);

%% Upstrokes (interpolated -20 mV upward crossings)
iu   = find(V(1:end-1) < v_thr & V(2:end) >= v_thr);
t_up = t(iu) + (v_thr - V(iu)) ./ (V(iu+1) - V(iu)) * dt_in;
nb   = numel(iu);
if nb < 4, error('Need at least 4 beats in %s (found %d).', infile, nb); end

%% MDP (min between successive upstrokes) and takeoff for each beat
% Beat k runs MDP(k) -> takeoff(k) -> upstroke(k+1) -> MDP(k+1)
iM = zeros(nb-1,1);
for k = 1:nb-1
    seg   = iu(k):iu(k+1);
    [~,j] = min(V(seg));
    iM(k) = seg(j);
end
nCyc = numel(iM) - 1;            % complete MDP-to-MDP cycles
tUp  = t_up(2:nCyc+1);           % upstroke that follows MDP(k)

% Takeoff (dV/dt criterion) -- used only for reporting
iTO = zeros(nCyc,1);
for k = 1:nCyc
    seg = iM(k):iu(k+1);
    j   = find(dVdt(seg) > dvdt_TO, 1, 'first');
    if isempty(j), j = 1; end
    iTO(k) = seg(j);
end
TOdelay = tUp - t(iTO);
fprintf('Takeoff -> upstroke delay: %.2f +/- %.2f ms\n', mean(TOdelay), std(TOdelay));

% Split each cycle at a fixed time BEFORE THE UPSTROKE, so the AP phase is
% upstroke-locked and the average keeps the full peak
pre = min(median(TOdelay), 0.8 * min(tUp - t(iM(1:nCyc))));
tS  = tUp - pre;                          % split time per beat
LA  = tS - t(iM(1:nCyc));                 % diastole (MDP -> split)
LB  = t(iM(2:nCyc+1)) - tS;               % AP (split -> next MDP), fixed
CL_ctrl = t(iM(2:nCyc+1)) - t(iM(1:nCyc));
fprintf('Control: %d cycles, CL = %.2f +/- %.2f ms (MDP-to-MDP)\n', ...
        nCyc, mean(CL_ctrl), std(CL_ctrl));
fprintf('  diastole %.2f ms, AP phase %.2f ms (starts %.2f ms before upstroke)\n', ...
        mean(LA), mean(LB), pre);

%% Ensemble-average each phase
% AP: absolute time from the split point (upstroke-locked), up to the shortest AP
tB = (0:dt_out:min(LB))';
B  = zeros(numel(tB), nCyc);
for k = 1:nCyc
    B(:,k) = interp1(t - tS(k), V, tB, 'linear');
end
pk_beats = max(B, [], 1);
B = mean(B, 2);
fprintf('Peak: beats %.2f +/- %.2f mV, template %.2f mV\n', mean(pk_beats), std(pk_beats), max(B));

% Diastole: normalized phase s in [0,1]
s  = linspace(0, 1, 2001)';
A  = zeros(numel(s), nCyc);
for k = 1:nCyc
    tk     = t(iM(k)) + s * LA(k);
    A(:,k) = interp1(t, V, tk, 'linear');
end
A = mean(A, 2);

% Make the wrap continuous: A must start where B ends (MDP)
% and end where B starts (split point)
A = A + (B(end) - A(1)) * (1 - s) + (B(1) - A(end)) * s;

%% Build and write templates
if ~exist(outdir, 'dir'), mkdir(outdir); end
tmpl = struct('CL', {}, 'V', {}, 'nA', {});
nB   = numel(B);
nNat = round(mean(LA) / dt_out);          % native diastole length (samples)
for c = 1:numel(CL_list)
    CL = CL_list(c);
    N  = round(CL / dt_out);
    nA = N - nB;
    if nA < 1
        error('CL = %g ms is shorter than the AP (%.1f ms).', CL, nB*dt_out);
    end

    if hold_mode && nA > nNat
        % Flat hold at the MDP, then the native (unstretched) diastolic ramp
        nHold = nA - nNat;
        sNat  = (0:nNat-1)' / nNat;       % excludes s = 1 (that sample is B(1))
        Vdia  = [repmat(A(1), nHold, 1); interp1(s, A, sNat, 'linear')];
        mode_str = sprintf('hold %.1f ms at MDP + native ramp %.1f ms', ...
                           nHold*dt_out, nNat*dt_out);
    else
        % Stretch (or compress) the diastolic ramp to fill the CL
        sA   = (0:nA-1)' / nA;            % excludes s = 1 (that sample is B(1))
        Vdia = interp1(s, A, sA, 'linear');
        mode_str = 'stretched ramp';
    end
    Vt = [Vdia; B];

    fn  = fullfile(outdir, sprintf('AP_%d.txt', CL));
    fid = fopen(fn, 'w');  fprintf(fid, '%.6f\n', Vt);  fclose(fid);
    fprintf('Wrote %s: %d samples (%s), MDP %.2f mV, peak %.2f mV\n', ...
            fn, numel(Vt), mode_str, min(Vt), max(Vt));

    tmpl(c).CL = CL;  tmpl(c).V = Vt;  tmpl(c).nA = nA;
end

%% Figure 1: control beats vs the template, aligned on the upstroke
% (for a long-CL hold template, the last ~native cycle is the part to compare)
figure('Color','w','Position',[100 100 900 380]);
subplot(1,2,1); hold on;
for k = 1:nCyc
    seg = iM(k):iM(k+1);
    plot(t(seg) - tUp(k), V(seg), 'Color', [0.75 0.75 0.75]);
end
[~, cNat] = min(abs(CL_list - mean(CL_ctrl)));
Vt = tmpl(cNat).V;  tt = ((0:numel(Vt)-1)' - tmpl(cNat).nA) * dt_out - pre;
plot(tt, Vt, 'k', 'LineWidth', 1.5);
xlim([-(mean(LA) + pre) - 10, mean(LB) - pre + 10]);
xlabel('Time from upstroke (ms)'); ylabel('V_m (mV)');
title(sprintf('Control beats (grey) vs template AP\\_%d (black)', CL_list(cNat)));
box off;

% Residual: template minus ensemble mean, on the AP phase
subplot(1,2,2); hold on;
Bk = zeros(nB, nCyc);
for k = 1:nCyc, Bk(:,k) = interp1(t - tS(k), V, tB, 'linear'); end
plot(tB - pre, Bk - B, 'Color', [0.75 0.75 0.75]);
yline(0, 'k');
xlabel('Time from upstroke (ms)'); ylabel('Beat - template (mV)');
title('Beat-to-beat deviation from template (AP phase)'); box off;

%% Figure 2: full template cycle(s), aligned on the upstroke
figure('Color','w','Position',[100 550 900 380]); hold on;
cols = lines(numel(tmpl));
for c = 1:numel(tmpl)
    Vt = tmpl(c).V;  tt = ((0:numel(Vt)-1)' - tmpl(c).nA) * dt_out - pre;
    plot(tt, Vt, 'Color', cols(c,:), 'LineWidth', 1.5, ...
         'DisplayName', sprintf('CL %d ms', tmpl(c).CL));
end
seg = iM(1):iM(2);
plot(t(seg) - tUp(1), V(seg), 'k--', 'LineWidth', 1, 'DisplayName', 'Control beat');
xlabel('Time from upstroke (ms)'); ylabel('V_m (mV)');
title('AP-clamp template(s), one cycle MDP to MDP');
legend('Location','northwest'); legend boxoff; box off;