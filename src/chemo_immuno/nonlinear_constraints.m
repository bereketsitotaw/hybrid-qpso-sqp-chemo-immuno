function [c, ceq] = nonlinear_constraints(z, p, N)
%NONLINEAR_CONSTRAINTS Trapezoidal collocation defects + boundary conditions.
%
% Defects follow Eq. (10) with Michaelis-Menten saturation in the RHS.

    idx = 1:(N + 1);
    x = z(idx);
    y = z(idx + (N + 1));
    z_vals = z(idx + 2 * (N + 1));
    u = z(idx + 3 * (N + 1));
    v = z(idx + 4 * (N + 1));
    T = z(end);
    h = T / N;

    fx = @(xx, yy, uu) p.xi * xx .* (1 - xx / p.x_inf) - p.theta * xx .* yy ...
        - p.alpha * xx .* (uu ./ (1 + uu));
    fy = @(xx, yy, zz, uu, vv) p.a * (1 - p.b * xx) .* yy .* zz + p.gamma - p.delta * yy ...
        - p.kappa * yy .* (uu ./ (1 + uu)) + p.nu * yy .* (vv ./ (1 + vv));
    fz = @(xx, zz, uu) p.sigma * xx + p.psi * xx .* (uu ./ (1 + uu)) - p.mu * zz;

    ceq_x = x(2:end) - x(1:end-1) - (h / 2) * ( ...
        fx(x(1:end-1), y(1:end-1), u(1:end-1)) + fx(x(2:end), y(2:end), u(2:end)));
    ceq_y = y(2:end) - y(1:end-1) - (h / 2) * ( ...
        fy(x(1:end-1), y(1:end-1), z_vals(1:end-1), u(1:end-1), v(1:end-1)) ...
        + fy(x(2:end), y(2:end), z_vals(2:end), u(2:end), v(2:end)));
    ceq_z = z_vals(2:end) - z_vals(1:end-1) - (h / 2) * ( ...
        fz(x(1:end-1), z_vals(1:end-1), u(1:end-1)) + fz(x(2:end), z_vals(2:end), u(2:end)));

    ceq_boundary = [x(1) - p.x0; y(1) - p.y0; z_vals(1) - p.z0; x(end) - p.xf];

    ceq = [ceq_x; ceq_y; ceq_z; ceq_boundary];
    c = [];
end
