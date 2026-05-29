%RUN_SENSITIVITY Reproduce Section 5.4 sensitivity / Pareto analysis.
clear; clc; close all;

addpath(fileparts(mfilename('fullpath')));
repo_root = init_project();

archive_script = fullfile(repo_root, 'archive', 'legacy', 'param_sensitivity_original.m');

if ~isfile(archive_script)
    error('Sensitivity script not found: %s', archive_script);
end

fprintf('Running archived sensitivity analysis (N=50 reduced grid)...\n');
fprintf('Source: %s\n\n', archive_script);
run(archive_script);

figures_src = fullfile(repo_root, 'results', 'figures');
if ~exist(figures_src, 'dir')
    mkdir(figures_src);
end
fprintf('\nNote: re-enable exportgraphics blocks in param_sensitivity_original.m to save EPS figures to results/figures/.\n');
