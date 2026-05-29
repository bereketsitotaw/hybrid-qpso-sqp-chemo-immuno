function cfg = apply_reproducibility_settings(cfg)
%APPLY_REPRODUCIBILITY_SETTINGS Pin RNG, ODE, and solver settings on a case cfg.

    repro = reproducibility_settings();

    rng(repro.rng_seed, repro.rng_generator);

    cfg.rng_seed = repro.rng_seed;
    cfg.rng_generator = repro.rng_generator;
    cfg.N = repro.N;

    cfg.qpso.swarmSize = repro.qpso_swarm_size;
    cfg.qpso.maxIter = repro.qpso_max_iter;
    cfg.qpso.rho = repro.qpso_rho;
    cfg.qpso.noise_amplitude = repro.qpso_noise_amplitude;
    cfg.qpso.ode_rel_tol = repro.ode_rel_tol;
    cfg.qpso.ode_abs_tol = repro.ode_abs_tol;

    cfg.sqp.algorithm = repro.sqp_algorithm;
    cfg.sqp.maxIterations = repro.sqp_max_iterations;
    cfg.sqp.maxFunctionEvaluations = repro.sqp_max_function_evaluations;
    cfg.sqp.constraintTolerance = repro.sqp_constraint_tolerance;
    cfg.sqp.optimalityTolerance = repro.sqp_optimality_tolerance;
    cfg.sqp.smooth_window = repro.sqp_smooth_window;

    cfg.repro = repro;
end
