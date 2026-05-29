%RUN_CASE_A Reproduce Section 5.1 (Case A: balanced objectives).
clear; clc; close all;

addpath(fileparts(mfilename('fullpath')));
repo_root = init_project();

cfg = apply_reproducibility_settings(case_a_balanced());

[solution, summary] = solve_chemo_immuno_ocp(cfg);
print_performance_summary(solution, summary);

results_dir = fullfile(repo_root, 'results');
save_solution(solution, summary, fullfile(results_dir, 'case_a_solution.mat'));
plot_chemo_immuno_results(solution, summary, fullfile(results_dir, 'figures'));
