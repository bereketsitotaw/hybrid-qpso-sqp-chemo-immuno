function repro = reproducibility_settings()
%REPRODUCIBILITY_SETTINGS Pinned environment for deterministic reproduction.
%
% Documents the author environment used to generate published figures.
% Table 2 in the paper lists RelTol/AbsTol = 1e-3 for QPSO warm-start;
% the pinned author settings below use tighter ODE tolerances for rollout
% stability during artifact reproduction.

    repro = struct();

    % --- Documented author environment (README / CITATION.cff) ---
    repro.matlab_version = 'R2024a';
    repro.optimization_toolbox = '9.5';  % update if your install differs
    repro.release_tag = 'v1.0.0-paper';

    % --- Random number generator (deterministic QPSO initialization) ---
    repro.rng_seed = 42;
    repro.rng_generator = 'twister';

    % --- QPSO rollout ODE tolerances (author pinned) ---
    repro.ode_rel_tol = 1e-8;
    repro.ode_abs_tol = 1e-10;

    % --- Paper Table 2 warm-start tolerances (documented alternative) ---
    repro.paper_table2_ode_rel_tol = 1e-3;
    repro.paper_table2_ode_abs_tol = 1e-3;

    % --- QPSO hyperparameters (Table 2) ---
    repro.qpso_swarm_size = 50;
    repro.qpso_max_iter = 100;
    repro.qpso_rho = 1e6;
    repro.qpso_noise_amplitude = 0.1;

    % --- SQP / collocation (Table 2) ---
    repro.sqp_algorithm = 'sqp';
    repro.sqp_max_iterations = 5000;
    repro.sqp_max_function_evaluations = 3e6;
    repro.sqp_constraint_tolerance = 1e-6;
    repro.sqp_optimality_tolerance = 1e-6;
    repro.sqp_smooth_window = 10;

    % --- Grid ---
    repro.N = 100;

    % --- Reference validation tolerances (relative to paper Table 3) ---
    repro.metric_abs_tol_T = 0.15;       % days
    repro.metric_abs_tol_exposure = 0.20; % normalized dose-time
    repro.metric_abs_tol_xf = 0.05;    % terminal tumor (target = 1.0)
end
