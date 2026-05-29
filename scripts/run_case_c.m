%RUN_CASE_C Reproduce Section 5.3 (Case C: minimum-time limit, MM kinetics).
clear; clc; close all;

addpath(fileparts(mfilename('fullpath')));
repo_root = init_project();

cfg = apply_reproducibility_settings(case_c_min_time());

[solution, summary] = solve_chemo_immuno_ocp(cfg);
print_performance_summary(solution, summary);

results_dir = fullfile(repo_root, 'results');
save_solution(solution, summary, fullfile(results_dir, 'case_c_solution.mat'));
plot_chemo_immuno_results(solution, summary, fullfile(results_dir, 'figures'));
