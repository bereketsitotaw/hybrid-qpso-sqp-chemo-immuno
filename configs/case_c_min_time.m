function cfg = case_c_min_time()
%CASE_C_MIN_TIME Section 5.3: minimum-time limit with MM kinetics (A=B=0, C>0).

    cfg = case_a_balanced();
    cfg.case_name = 'Case C: Minimum-time limit';
    cfg.p.l1 = 0;
    cfg.p.l2 = 0;
    cfg.p.l3 = 0;
    cfg.p.A = 0.0;
    cfg.p.B = 0.0;
    cfg.p.C = 1.0;
    cfg.T_guess = 3.0;
    cfg.sqp.clamp_bang_start = false;
end
