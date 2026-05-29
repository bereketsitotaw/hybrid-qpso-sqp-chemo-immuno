%REPRODUCE_FIGURE_5 Bound sweep for Case A (Paper Figure 5).
%
% Status: uses archived sensitivity workflow (partial). Full umax/vmax sweep
% driver is not yet modularized; this script documents the gap and runs the
% closest available archived analysis.
fprintf(['Figure 5 (bound sweep) requires a dedicated sweep driver.\n', ...
    'Running archived sensitivity script as partial proxy (weight sweeps, not bounds).\n', ...
    'See docs/FIGURE_MAP.md for status.\n\n']);
addpath(fileparts(mfilename('fullpath')));
repo_root = init_project();
run(fullfile(repo_root, 'scripts', 'run_sensitivity.m'));
