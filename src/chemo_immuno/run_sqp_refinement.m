function [z_opt, J_opt, lambda_struct, fmincon_time, sqp] = run_sqp_refinement(p, N, gBest, T_guess, sqp)
%RUN_SQP_REFINEMENT Stage-2 direct collocation with fmincon SQP (Section 4).

    u_guess = gBest(1:N + 1);
    v_guess = gBest(N + 2:end);

    if sqp.smooth_controls
        u_guess = smoothdata(u_guess, 'gaussian', sqp.smooth_window);
        v_guess = smoothdata(v_guess, 'gaussian', sqp.smooth_window);
    end

    if sqp.clamp_bang_start && u_guess(1) > 0.8 * p.umax
        u_guess(1) = p.umax;
    end

    x_guess = linspace(p.x0, p.xf, N + 1);
    y_guess = p.y0 * ones(1, N + 1);
    z_guess = p.z0 * ones(1, N + 1);
    z0 = [x_guess, y_guess, z_guess, u_guess, v_guess, T_guess]';

    lb_x = zeros(1, N + 1);
    ub_x = p.x_inf * ones(1, N + 1);
    lb_y = zeros(1, N + 1);
    ub_y = inf(1, N + 1);
    lb_z = zeros(1, N + 1);
    ub_z = inf(1, N + 1);
    lb_u = zeros(1, N + 1);
    ub_u = p.umax * ones(1, N + 1);
    lb_v = zeros(1, N + 1);
    ub_v = p.vmax * ones(1, N + 1);
    lb_T = sqp.T_min;
    ub_T = sqp.T_max;

    lb = [lb_x, lb_y, lb_z, lb_u, lb_v, lb_T]';
    ub = [ub_x, ub_y, ub_z, ub_u, ub_v, ub_T]';

    options = optimoptions('fmincon', ...
        'Algorithm', sqp.algorithm, ...
        'Display', sqp.display, ...
        'MaxIterations', sqp.maxIterations, ...
        'MaxFunctionEvaluations', sqp.maxFunctionEvaluations, ...
        'ConstraintTolerance', sqp.constraintTolerance, ...
        'OptimalityTolerance', sqp.optimalityTolerance);

    tic;
    [z_opt, J_opt, ~, ~, lambda_struct] = fmincon( ...
        @(z) objective_func(z, p, N), z0, ...
        [], [], [], [], lb, ub, @(z) nonlinear_constraints(z, p, N), options);
    fmincon_time = toc;
end
