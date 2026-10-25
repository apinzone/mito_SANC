function compare_models(backup_mat, updated_mat)
%COMPARE_MODELS  Load two analyze_model() .mat result files and plot a
%   Vm/Cai overlay comparison. Run this after analyze_model() has been run
%   once in each model's own directory.
%
%   compare_models('Backup_results.mat', 'Updated_results.mat')

b = load(backup_mat);  backup  = b.results;
u = load(updated_mat); updated = u.results;

figure('Name', 'Vm and Cai Comparison', 'Position', [100 100 900 700], 'Color', 'w');

subplot(2,1,1);
plot(backup.t - backup.t(1), backup.vm, 'LineWidth', 1.0); hold on;
plot(updated.t - updated.t(1), updated.vm, 'LineWidth', 1.0);
ylabel('V_m (mV)');
title('Whole-Cell V_m and Ca_i Comparison (last analysis window)');
legend({backup.label, updated.label}, 'Location', 'best');

subplot(2,1,2);
plot(backup.t - backup.t(1), backup.cai, 'LineWidth', 1.0); hold on;
plot(updated.t - updated.t(1), updated.cai, 'LineWidth', 1.0);
ylabel('Ca_i (\muM)');
xlabel('Time since window start (ms)');
legend({backup.label, updated.label}, 'Location', 'best');

saveas(gcf, 'vm_cai_comparison.png');
fprintf('\nSaved comparison plot: vm_cai_comparison.png\n');

% --- inak / Iup comparison (both models compute these, ATP-dependent or not) ---
if isfield(backup.currents, 'inak') && isfield(updated.currents, 'inak') ...
        && isfield(backup.currents, 'Iup_avg') && isfield(updated.currents, 'Iup_avg')

    tb = backup.currents.time;  tu = updated.currents.time;

    figure('Name', 'inak and Iup Comparison', 'Position', [150 150 900 700], 'Color', 'w');

    subplot(2,1,1);
    plot(tb - tb(1), backup.currents.inak, 'LineWidth', 1.0); hold on;
    plot(tu - tu(1), updated.currents.inak, 'LineWidth', 1.0);
    ylabel('I_{NaK} (pA)');
    title('I_{NaK} and I_{up} Comparison (last analysis window)');
    legend({backup.label, updated.label}, 'Location', 'best');

    subplot(2,1,2);
    plot(tb - tb(1), backup.currents.Iup_avg, 'LineWidth', 1.0); hold on;
    plot(tu - tu(1), updated.currents.Iup_avg, 'LineWidth', 1.0);
    ylabel('I_{up,avg} (\muM/ms)');
    xlabel('Time since window start (ms)');
    legend({backup.label, updated.label}, 'Location', 'best');

    saveas(gcf, 'inak_iup_comparison.png');
    fprintf('Saved comparison plot: inak_iup_comparison.png\n');
else
    fprintf('\nSkipping inak/Iup comparison plot: one or both models are missing a currents trace.\n');
    fprintf('(Backup needs inak_iup_trace.txt passed as currents_path in analyze_model().)\n');
end

% --- Quick side-by-side biomarker printout ---
fprintf('\n=== Side-by-side summary ===\n');
print_pair('BCL (ms)', backup.BCL_vm, updated.BCL_vm);
fn = fieldnames(backup.cat_markers);
for i = 1:numel(fn)
    if isfield(updated.cat_markers, fn{i})
        print_pair(fn{i}, backup.cat_markers.(fn{i}), updated.cat_markers.(fn{i}));
    end
end

end

function print_pair(name, a, b)
a = a(~isnan(a)); b = b(~isnan(b));
fprintf('  %-16s  %s: %.4f +/- %.4f    %s: %.4f +/- %.4f\n', ...
    name, 'backup', mean(a), std(a), 'updated', mean(b), std(b));
end