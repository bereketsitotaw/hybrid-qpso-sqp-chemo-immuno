function pmp = compute_pmp_diagnostics(x_opt, y_opt, z_opt, u_opt, v_opt, T_opt, N, p, lambda_struct)
%COMPUTE_PMP_DIAGNOSTICS PMP/KKT verification from collocation multipliers.

    time_opt = linspace(0, T_opt, N + 1);

    lambda_multipliers = lambda_struct.eqnonlin;
    h = T_opt / N;
    lambda1 = -lambda_multipliers(1:N) / h;
    lambda2 = -lambda_multipliers(N + 1:2 * N) / h;
    lambda3 = -lambda_multipliers(2 * N + 1:3 * N) / h;

    lambda1_opt = [lambda1; lambda1(end)];
    lambda2_opt = [lambda2; lambda2(end)];
    lambda3_opt = [lambda3; lambda3(end)];

    switching_term_u = -lambda1_opt .* p.alpha .* x_opt ...
        - lambda2_opt .* p.kappa .* y_opt ...
        + lambda3_opt .* p.psi .* x_opt;
    switching_term_v = lambda2_opt .* p.nu .* y_opt;

    dH_du = p.A + switching_term_u ./ (1 + u_opt).^2;
    dH_dv = p.B + switching_term_v ./ (1 + v_opt).^2;

    d2H_du2 = -2 * switching_term_u ./ (1 + u_opt).^3;
    d2H_dv2 = -2 * switching_term_v ./ (1 + v_opt).^3;

    dt = T_opt / N;
    dH_du_dt = gradient(dH_du, dt);
    dH_dv_dt = gradient(dH_dv, dt);

    H_u_contribution = p.A * u_opt + switching_term_u .* (u_opt ./ (1 + u_opt));
    H_v_contribution = p.B * v_opt + switching_term_v .* (v_opt ./ (1 + v_opt));
    H_control_contribution = H_u_contribution + H_v_contribution;

    H_state_lagrangian = p.l1 * x_opt + p.l2 * y_opt + p.l3 * z_opt + p.C;
    H_state_dynamics = lambda1_opt .* (p.xi * x_opt .* (1 - x_opt / p.x_inf) - p.theta * x_opt .* y_opt) ...
        + lambda2_opt .* (p.a * (1 - p.b * x_opt) .* y_opt .* z_opt + p.gamma - p.delta * y_opt) ...
        + lambda3_opt .* (p.sigma * x_opt - p.mu * z_opt);
    H_state_contribution = H_state_lagrangian + H_state_dynamics;
    H_total = H_state_contribution + H_control_contribution;

    pmp = struct();
    pmp.time = time_opt;
    pmp.lambda1 = lambda1_opt;
    pmp.lambda2 = lambda2_opt;
    pmp.lambda3 = lambda3_opt;
    pmp.switching_term_u = switching_term_u;
    pmp.switching_term_v = switching_term_v;
    pmp.dH_du = dH_du;
    pmp.dH_dv = dH_dv;
    pmp.d2H_du2 = d2H_du2;
    pmp.d2H_dv2 = d2H_dv2;
    pmp.dH_du_dt = dH_du_dt;
    pmp.dH_dv_dt = dH_dv_dt;
    pmp.H_state_contribution = H_state_contribution;
    pmp.H_control_contribution = H_control_contribution;
    pmp.H_total = H_total;
end
