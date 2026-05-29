function qpso = default_qpso_settings()
%DEFAULT_QPSO_SETTINGS QPSO hyperparameters from Table 2 (paper).

    qpso = struct();
    qpso.swarmSize = 50;
    qpso.maxIter = 100;
    qpso.rho = 1e6;
    qpso.noise_amplitude = 0.1;
    qpso.useParfor = false;
    qpso.verbose = true;
    qpso.ode_rel_tol = 1e-3;   % paper Table 2 default
    qpso.ode_abs_tol = 1e-3;   % overridden by apply_reproducibility_settings
end
