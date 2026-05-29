%REPRODUCE_FIGURE_3 Case B state/control + PMP (Paper Figure 3).
addpath(fileparts(mfilename('fullpath')));
repo_root = init_project();
cfg = apply_reproducibility_settings(case_b_asymmetric());
[solution, summary] = solve_chemo_immuno_ocp(cfg);
out_dir = fullfile(repo_root, 'results', 'figures');
export_publication_figures(solution, 'fig3_caseB', out_dir, {'png', 'eps'});
export_publication_figures(solution, 'fig3_pmp_caseB', out_dir, {'png', 'eps'});
save_solution(solution, summary, fullfile(repo_root, 'results', 'case_b_solution.mat'));
fprintf('Figure 3 outputs written to %s\n', out_dir);
