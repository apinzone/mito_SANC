function S = poster_style()
%POSTER_STYLE  One place for every poster-figure setting, plus two helpers.
%
%   S = poster_style();
%   fig = S.newfig(16, 10);            % 16 x 10 cm figure (= the poster slot size)
%   t = tiledlayout(fig, 2, 2, 'TileSpacing','compact', 'Padding','compact');
%   ax = nexttile(t); plot(ax, x, y, 'Color', S.c.orig, 'LineWidth', S.lw);
%   ... (all plotting) ...
%   S.finish(fig);                     % applies fonts/ticks to every axis, turns on edit mode
%   S.export(fig, 'fig2_control');     % writes .pdf (vector), .png (600 dpi), .fig
%
% Manual edits: after S.finish the figure is in plot-edit mode (click to select,
% double-click for properties). Resize the window if it helps you work - S.export
% always snaps the figure back to its poster size (stored in fig.UserData) first,
% so the exported files keep the right size and font scale.
% Re-export after editing:  S.export(gcf, 'name')
% Edit a saved figure later: openfig('name.fig'); S = poster_style(); ... S.export(gcf,'name')
%
% Rules that make figures look clean on the poster:
%   1. Set the FINAL physical size here (cm), equal to the box on the poster.
%      Never stretch the figure in MATLAB, PowerPoint or Illustrator afterwards:
%      stretching is what distorts fonts and line widths.
%   2. Fonts in points at that final size. Every figure uses the same sizes.
%   3. Same colour for the same model/condition in every figure on the poster.
%   4. No titles inside figures - the poster panel heading does that job.
%   5. Export vector PDF and place it at 100% scale.

% ---- fonts / lines (pt). For a 48x36 in poster viewed at ~1-2 m:
S.font      = 'Arial';
S.fs_label  = 20;        % axis labels
S.fs_tick   = 18;        % tick labels
S.fs_panel  = 24;        % panel letters A, B, C ...
S.lw        = 2.0;       % data lines
S.lw_axis   = 1.2;       % axis lines
S.edit_mode = true;      % open figures in plot-edit mode after S.finish

% ---- one colour per model / condition, used everywhere
S.c.orig    = [0.35 0.35 0.35];   % original model (no ATP)       - grey
S.c.mito    = [0.00 0.45 0.74];   % mito model, original params  - blue
S.c.tuned   = [0.85 0.33 0.10];   % mito model, tuned 6/6        - orange
S.c.song    = [0.47 0.67 0.19];   % Song et al.                  - green
S.c.mcu     = [0.30 0.45 0.80];   % MCU block (shade lighter for 50%)
S.c.nclx    = [0.90 0.55 0.15];   % NCLX block

S.newfig = @newfig;
S.finish = @(fig) finish(fig, S);
S.export = @export;
end

% ------------------------------------------------------------------------
function fig = newfig(w_cm, h_cm)
fig = figure('Color','w', 'Units','centimeters', 'Position',[2 2 w_cm h_cm], ...
             'PaperUnits','centimeters', 'PaperSize',[w_cm h_cm], ...
             'PaperPosition',[0 0 w_cm h_cm]);
fig.UserData.poster_size_cm = [w_cm h_cm];       % remembered by export, saved in the .fig
end

function finish(fig, S)
ax = findall(fig, 'Type','axes');
set(ax, 'FontName',S.font, 'FontSize',S.fs_tick, 'LineWidth',S.lw_axis, ...
        'Box','off', 'TickDir','out', 'TickLength',[0.02 0.02], 'Color','none');
for a = ax'
    a.XLabel.FontSize = S.fs_label;  a.YLabel.FontSize = S.fs_label;
    a.Title.String = '';                                   % rule 4
end
lg = findall(fig, 'Type','legend');
set(lg, 'FontName',S.font, 'FontSize',S.fs_tick, 'Box','off');
if S.edit_mode, plotedit(fig, 'on'); end
end

function export(fig, name)
% Snap back to the poster size before writing, whatever happened in the GUI.
was_editing = plotedit(fig, 'isactive');
if isstruct(fig.UserData) && isfield(fig.UserData, 'poster_size_cm')
    sz = fig.UserData.poster_size_cm;
    plotedit(fig, 'off');                                  % deselect so handles don't print
    fig.Units = 'centimeters';  fig.Position(3:4) = sz;
    fig.PaperUnits = 'centimeters';  fig.PaperSize = sz;  fig.PaperPosition = [0 0 sz];
    drawnow;
end
exportgraphics(fig, [name '.pdf'], 'ContentType','vector');   % for the poster
exportgraphics(fig, [name '.png'], 'Resolution',600);         % for slides / email
savefig(fig, [name '.fig']);                                  % for later edits
if was_editing, plotedit(fig, 'on'); end                      % back to edit mode
fprintf('Exported %s (.pdf/.png/.fig)\n', name);
end
