function cfg = case_b_asymmetric()
%CASE_B_ASYMMETRIC Section 5.2: A=1, B=0, C=1 (asymmetric regularization).

    cfg = case_a_balanced();
    cfg.case_name = 'Case B: Asymmetric regularization';
    cfg.p.B = 0.0;
end
