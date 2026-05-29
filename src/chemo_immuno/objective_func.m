function J = objective_func(z, p, N)
%OBJECTIVE_FUNC Trapezoidal quadrature of the Bolza running cost (Section 2.2).

    idx = 1:(N + 1);
    x = z(idx);
    y = z(idx + (N + 1));
    z_vals = z(idx + 2 * (N + 1));
    u = z(idx + 3 * (N + 1));
    v = z(idx + 4 * (N + 1));
    T = z(end);

    time = linspace(0, T, N + 1);
    integrand = p.l1 * x' + p.l2 * y' + p.l3 * z_vals' + p.A * u' + p.B * v' + p.C;
    J = trapz(time, integrand);
end
