function cost = qpso_cost_func(u, v, T, p, N, rho, qpso)
%QPSO_COST_FUNC Augmented objective for QPSO warm start (Eq. 9).

    if nargin < 6 || isempty(rho)
        rho = 1e6;
    end
    if nargin < 7
        qpso = struct();
    end

    time = linspace(0, T, N + 1);
    w0 = [p.x0; p.y0; p.z0];

    rel_tol = 1e-3;
    abs_tol = 1e-3;
    if isfield(qpso, 'ode_rel_tol') && ~isempty(qpso.ode_rel_tol)
        rel_tol = qpso.ode_rel_tol;
    end
    if isfield(qpso, 'ode_abs_tol') && ~isempty(qpso.ode_abs_tol)
        abs_tol = qpso.ode_abs_tol;
    end

    options = odeset('RelTol', rel_tol, 'AbsTol', abs_tol);
    [~, w_traj] = ode45(@(t, w) dynamics(t, w, u, v, p, time), time, w0, options);

    x = w_traj(:, 1)';
    y = w_traj(:, 2)';
    z_vals = w_traj(:, 3)';

    integrand = p.l1 * x + p.l2 * y + p.l3 * z_vals + p.A * u + p.B * v + p.C;
    cost = trapz(time, integrand);
    cost = cost + rho * (x(end) - p.xf)^2;
end
