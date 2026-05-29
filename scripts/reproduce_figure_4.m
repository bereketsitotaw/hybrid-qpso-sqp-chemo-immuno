%REPRODUCE_FIGURE_4 Case C minimum-time limit (Paper Figure 4).
addpath(fileparts(mfilename('fullpath')));
repo_root = init_project();
cfg = apply_reproducibility_settings(case_c_min_time());
[solution, summary] = solve_chemo_immuno_ocp(cfg);
out_dir = fullfile(repo_root, 'results', 'figures');
export_publication_figures(solution, 'fig4_caseC', out_dir, {'png', 'eps'});
export_publication_figures(solution, 'fig4_pmp_caseC', out_dir, {'png', 'eps'});
save_solution(solution, summary, fullfile(repo_root, 'results', 'case_c_solution.mat'));
fprintf('Figure 4 outputs written to %s\n', out_dir);
