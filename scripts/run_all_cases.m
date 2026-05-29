%RUN_ALL_CASES Execute Cases A, B, and C with pinned reproducibility settings.
addpath(fileparts(mfilename('fullpath')));
init_project();
reproduce_paper_pipeline('export_figures', true, 'run_sensitivity', false);
