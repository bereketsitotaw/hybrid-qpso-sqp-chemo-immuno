function dwdt = dynamics(t, w, u_vec, v_vec, p, time_vec)
%DYNAMICS Michaelis-Menten chemo-immunotherapy state equations (Eqs. 1-3).
%
% States: w = [x; y; z]  (tumor, immune effector, antigen)
% Controls are linearly interpolated from discrete node values.

    x = w(1);
    y = w(2);
    z = w(3);

    u = interp1(time_vec, u_vec, t, 'linear', 'extrap');
    v = interp1(time_vec, v_vec, t, 'linear', 'extrap');

    u_sat = u / (1 + u);
    v_sat = v / (1 + v);

    dxdt = p.xi * x * (1 - x / p.x_inf) - p.theta * x * y - p.alpha * x * u_sat;
    dydt = p.a * (1 - p.b * x) * y * z + p.gamma - p.delta * y ...
           - p.kappa * y * u_sat + p.nu * y * v_sat;
    dzdt = p.sigma * x + p.psi * x * u_sat - p.mu * z;

    dwdt = [dxdt; dydt; dzdt];
end
