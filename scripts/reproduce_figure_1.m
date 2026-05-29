%REPRODUCE_FIGURE_1 Case A state/control figure (Paper Figure 1).
addpath(fileparts(mfilename('fullpath')));
repo_root = init_project();
cfg = apply_reproducibility_settings(case_a_balanced());
[solution, summary] = solve_chemo_immuno_ocp(cfg);
out_dir = fullfile(repo_root, 'results', 'figures');
export_publication_figures(solution, 'fig1_caseA', out_dir, {'png', 'eps'});
save_solution(solution, summary, fullfile(repo_root, 'results', 'case_a_solution.mat'));
fprintf('Figure 1 outputs written to %s\n', out_dir);
