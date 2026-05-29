%VALIDATE_REPO Lightweight integrity checks (run in MATLAB).
addpath(fileparts(mfilename('fullpath')));
repo_root = init_project();

required = {
    fullfile(repo_root, 'paper', 'published_paper.pdf')
    fullfile(repo_root, 'src', 'chemo_immuno', 'solve_chemo_immuno_ocp.m')
    fullfile(repo_root, 'configs', 'case_a_balanced.m')
    fullfile(repo_root, 'configs', 'reproducibility_settings.m')
    fullfile(repo_root, 'scripts', 'reproduce_paper.m')
    fullfile(repo_root, 'docs', 'FIGURE_MAP.md')
    fullfile(repo_root, 'CITATION.cff')
    fullfile(repo_root, 'archive', 'legacy', 'qpso_optimal_original.m')
    };

for i = 1:numel(required)
    assert(isfile(required{i}), 'Missing required file: %s', required{i});
end

fprintf('Required files: OK\n');

cfg = apply_reproducibility_settings(case_a_balanced());
assert(cfg.N == 100, 'Expected N=100 for paper grid.');
assert(cfg.rng_seed == 42, 'Expected pinned rng seed 42.');
assert(cfg.qpso.ode_rel_tol == 1e-8, 'Expected pinned ODE RelTol 1e-8.');

fprintf('Reproducibility config: OK\n');

test_dynamics();
fprintf('\nRepository validation complete.\n');
fprintf('Next: run(''scripts/reproduce_paper.m'') for full reproduction.\n');
