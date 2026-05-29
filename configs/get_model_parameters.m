function p = get_model_parameters()
%GET_MODEL_PARAMETERS Biological parameters from Table 1 (paper).

    p = struct();
    p.xi = 0.5618;
    p.x_inf = 780;
    p.theta = 1.0;
    p.alpha = 1.0;
    p.a = 0.00726;
    p.b = 0.00264;
    p.gamma = 0.1181;
    p.delta = 0.3745;
    p.kappa = 2.0;
    p.nu = 1.0;
    p.sigma = 0.1;
    p.psi = 2.0;
    p.mu = 0.45;

    p.x0 = 650;
    p.y0 = 0.1;
    p.z0 = 400;
    p.xf = 1.0;

    p.umax = 1.0;
    p.vmax = 1.0;
end
