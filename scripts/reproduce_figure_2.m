%REPRODUCE_FIGURE_2 Case A PMP diagnostics (Paper Figure 2).
addpath(fileparts(mfilename('fullpath')));
repo_root = init_project();
cfg = apply_reproducibility_settings(case_a_balanced());
[solution, summary] = solve_chemo_immuno_ocp(cfg);
out_dir = fullfile(repo_root, 'results', 'figures');
export_publication_figures(solution, 'fig2_pmp_caseA', out_dir, {'png', 'eps'});
fprintf('Figure 2 outputs written to %s\n', out_dir);
