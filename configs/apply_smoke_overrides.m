function cfg = apply_smoke_overrides(cfg)
%APPLY_SMOKE_OVERRIDES Tiny non-paper QPSO/SQP settings for plumbing checks.

    repro = reproducibility_settings();

    cfg.N = repro.smoke_grid_N;
    cfg.qpso.swarmSize = repro.smoke_qpso_swarm_size;
    cfg.qpso.maxIter = repro.smoke_qpso_max_iter;
    cfg.qpso.useParfor = repro.smoke_use_parfor;
    cfg.qpso.verbose = false;

    cfg.sqp.maxIterations = repro.smoke_sqp_max_iterations;
    cfg.sqp.maxFunctionEvaluations = repro.smoke_sqp_max_function_evaluations;
    cfg.sqp.display = repro.smoke_sqp_display;

    cfg.run_pmp = repro.smoke_enable_pmp;
    cfg.smoke = true;
    cfg.case_name = [cfg.case_name ' [SMOKE non-paper]'];
end
