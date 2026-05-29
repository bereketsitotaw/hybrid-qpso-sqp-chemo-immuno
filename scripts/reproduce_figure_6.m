%REPRODUCE_FIGURE_6 Parameter sensitivity tornado (Paper Figure 6).
%
% Status: partial — archived param_sensitivity includes weight/Pareto panels
% but not the full ±20%% biological parameter re-optimization tornado.
fprintf(['Figure 6 (parameter tornado) is partially covered by param_sensitivity_original.m.\n', ...
    'See docs/FIGURE_MAP.md. Running archived sensitivity script...\n\n']);
addpath(fileparts(mfilename('fullpath')));
repo_root = init_project();
run(fullfile(repo_root, 'archive', 'legacy', 'param_sensitivity_original.m'));
