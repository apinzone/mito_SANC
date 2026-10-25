%% make_fig3_linescan.m
% Fig. 3 with spatial detail: 2 x 3 panels, all on the same 3 aligned beats.
%
%     Vm        I_NaK      mito Ca linescan - peripheral row
%     [Ca]i     J_up       mito Ca linescan - central row
%
% Traces: original model (grey) vs ATP-dependent model (blue dashed), as in Fig. 3.
% Linescans: ATP-dependent model only, from the SAME run as its traces, so every
% stripe lines up with a beat in the traces.
%
% snap_dir must contain (all from the patched SAN_pace.cpp run):
%   ci.txt, ikatp_inak_serca_trace.txt, mito_map.txt, mito_ca_grid.txt
% orig_dir: ci.txt + inak_iup_trace.txt (original model, as before)
% Needs poster_style.m. Output: fig3_linescan.pdf/.png/.fig

clear; close all;
S = poster_style();

%% ===================== SETTINGS =====================
base     = pwd;
orig_dir = fullfile(base, 'orig_dir');
snap_dir = fullfile(base, 'snap_dir');    % <-- the new 60 s run with snapshot output
fig_size = [46 20.5];                     % cm
nBeats   = 3;
pre      = 45;                            % ms before the first aligned upstroke
post     = 35;                            % ms trimmed from the end
row_per  = [0 1];                         % [my mz] peripheral mito row
row_cen  = [1 1];                         % [my mz] central mito row
dx_um    = 1.84;                          % longitudinal CRU spacing (um)
smooth_f = 4;                             % display upsampling along position (1 = none)
fname    = 'fig3_linescan';
% =====================================================

runs = struct('dir',  {orig_dir,   snap_dir}, ...
              'name', {'Original', 'ATP-dependent'}, ...
              'color',{S.c.orig,   S.c.mito}, ...
              'style',{'-',        '--'});
keys = {'Vm','Cai','INaK','Iup'};

fig = S.newfig(fig_size(1), fig_size(2));
tl  = tiledlayout(fig, 2, 3, 'TileSpacing','compact', 'Padding','compact');
ax  = gobjects(6,1);
tile = [1 4 2 5 3 6];                     % column-major: Vm, Cai | INaK, Iup | linescans
for i = 1:6, ax(i) = nexttile(tl, tile(i)); hold(ax(i),'on'); end

%% ---- traces (same alignment rule as Fig. 3) ----
CLref = [];  t0_snap = [];
for k = 1:numel(runs)
    D  = load_run(runs(k).dir);
    tu = upstrokes(D.t, D.vm);
    if isempty(CLref), CLref = median(diff(tu(end-nBeats:end))); end
    t0 = tu(end-nBeats);
    if strcmp(runs(k).dir, snap_dir), t0_snap = t0; end
    for i = 1:4
        [tt, yy, lab] = getvar(D, keys{i});
        if isempty(tt), continue; end
        w = tt >= t0 - pre & tt <= t0 + nBeats*CLref - post;
        plot(ax(i), tt(w) - t0, yy(w), runs(k).style, 'Color',runs(k).color, ...
             'LineWidth',S.lw, 'DisplayName',runs(k).name);
        ylabel(ax(i), lab);
    end
end
xl = [-pre, nBeats*CLref - post];

%% ---- linescans (ATP-dependent run, same t0 as its traces) ----
M  = readtable(fullfile(snap_dir,'mito_map.txt'), 'FileType','text', 'Delimiter','\t');
CA = readmatrix(fullfile(snap_dir,'mito_ca_grid.txt'), 'FileType','text', 'Delimiter','\t');
tg = CA(:,1) - t0_snap;  ca = CA(:,2:end)*1000;          % nM
w  = tg >= xl(1) & tg <= xl(2);
if ~any(w), error('Snapshot window does not cover the aligned beats - increase snap_window_ms.'); end
tg = tg(w);  ca = ca(w,:);
dca = ca - min(ca, [], 1);                                % delta mito Ca (each mito vs its own min)

rows = {row_per, row_cen};
for r = 1:2
    ii = find(M.my==rows{r}(1) & M.mz==rows{r}(2));       % rows of M = columns of ca
    [~, o] = sort(M.mx(ii));  ii = ii(o);
    img = dca(:, ii)';                                    % position x time
    pos = (0:numel(ii)-1) * dx_um;
    if smooth_f > 1                                       % display-only upsampling along position
        posf = linspace(pos(1), pos(end), (numel(pos)-1)*smooth_f + 1);
        img  = interp1(pos, img, posf, 'linear');  pos = posf;
    end
    rows{r} = struct('img', img, 'pos', pos, 'ii', ii);
end
cl = [0 max([rows{1}.img(:); rows{2}.img(:)])];          % shared colour scale
ttl = {'Peripheral mitochondria','Central mitochondria'};
for r = 1:2
    a = ax(4+r);
    imagesc(a, tg, rows{r}.pos, rows{r}.img);  axis(a, 'xy');  caxis(a, cl);  colormap(a, parula(256));
    ylabel(a, 'Position (\mum)');
    text(a, 0.5, 1.04, ttl{r}, 'Units','normalized', 'HorizontalAlignment','center', ...
         'FontSize',S.fs_label, 'FontName',S.font, 'FontWeight','bold');
end
cb = colorbar(ax(6));  cb.Label.String = '\Delta[Ca^{2+}]_{mito} (nM)';
fprintf('Mean beat swing: peripheral row %.1f nM | central row (ends excluded) %.1f nM\n', ...
    mean(max(dca(:, rows{1}.ii))), mean(max(dca(:, rows{2}.ii(3:end-2)))));

%% ---- axes, labels, export ----
for i = 1:6
    if any(i == [2 4 6]), xlabel(ax(i), 'Time (ms)'); else, set(ax(i), 'XTickLabel', []); end
    text(ax(i), -0.15, 1.04, char('A'+i-1), 'Units','normalized', ...
         'FontWeight','bold', 'FontSize',S.fs_panel, 'FontName',S.font);
end
linkaxes(ax, 'x');  xlim(ax(1), xl);
legend(ax(1), 'Location','southwest');
S.finish(fig);
set(ax(5:6), 'Layer','top');
cb.FontSize = S.fs_tick;  cb.FontName = S.font;  cb.Label.FontSize = S.fs_label;
S.export(fig, fname);

%% =====================================================================
function D = load_run(d)
C = readmatrix(fullfile(d,'ci.txt'), 'FileType','text', 'Delimiter',' ');
D.t = C(:,1);  D.vm = C(:,2);  D.cai = C(:,3);
D.T = [];
if isfile(fullfile(d,'ikatp_inak_serca_trace.txt'))
    D.T = readtable(fullfile(d,'ikatp_inak_serca_trace.txt'), 'FileType','text', 'Delimiter','\t');
elseif isfile(fullfile(d,'inak_iup_trace.txt'))
    D.T = readtable(fullfile(d,'inak_iup_trace.txt'), 'FileType','text', 'Delimiter','\t');
end
end

function tu = upstrokes(t, vm)
iu = find(vm(1:end-1) < -20 & vm(2:end) >= -20);
tu = t(iu) + (-20 - vm(iu)) ./ (vm(iu+1) - vm(iu)) .* (t(iu+1) - t(iu));
end

function [tt, yy, lab] = getvar(D, key)
tt = []; yy = []; lab = '';
switch key
    case 'Vm',   tt = D.t; yy = D.vm;        lab = 'V_m (mV)';
    case 'Cai',  tt = D.t; yy = 1000*D.cai;  lab = '[Ca^{2+}]_i (nM)';
    case 'INaK', if ~isempty(D.T), tt = D.T.time; yy = D.T.inak/0.025; lab = 'I_{NaK} (pA/pF)'; end
    case 'Iup',  if ~isempty(D.T), tt = D.T.time; yy = D.T.Iup_avg;    lab = 'J_{up} (\muM/ms)'; end
end
end