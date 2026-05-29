function cfg = case_a_balanced()
%CASE_A_BALANCED Section 5.1: (l,A,B,C) = (0,0,0; 1, 0.3, 1).

    cfg = struct();
    cfg.case_name = 'Case A: Balanced objectives';
    cfg.p = get_model_parameters();
    cfg.p.l1 = 0;
    cfg.p.l2 = 0;
    cfg.p.l3 = 0;
    cfg.p.A = 1.0;
    cfg.p.B = 0.3;
    cfg.p.C = 1.0;

    cfg.N = 100;
    cfg.T_guess = 4.0;
    cfg.qpso = default_qpso_settings();
    cfg.sqp = default_sqp_settings();
    cfg.run_pmp = true;
    cfg.rng_seed = 42;
end
