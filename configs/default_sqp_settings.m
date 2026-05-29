function sqp = default_sqp_settings()
%DEFAULT_SQP_SETTINGS SQP/collocation settings from Table 2 (paper).

    sqp = struct();
    sqp.algorithm = 'sqp';
    sqp.display = 'iter';
    sqp.maxIterations = 5000;
    sqp.maxFunctionEvaluations = 3e6;
    sqp.constraintTolerance = 1e-6;
    sqp.optimalityTolerance = 1e-6;
    sqp.smooth_controls = true;
    sqp.smooth_window = 10;
    sqp.clamp_bang_start = true;
    sqp.T_min = 0.1;
    sqp.T_max = 20;
end
