function test_case_weights()
%TEST_CASE_WEIGHTS Section 5 case (A,B,C) match paper captions.

    a = case_a_balanced();
    assert(abs(a.p.A - 1.0) < 1e-12 && abs(a.p.B - 0.3) < 1e-12 && abs(a.p.C - 1.0) < 1e-12);
    assert(a.p.l1 == 0 && a.p.l2 == 0 && a.p.l3 == 0);

    b = case_b_asymmetric();
    assert(abs(b.p.A - 1.0) < 1e-12 && abs(b.p.B - 0.0) < 1e-12 && abs(b.p.C - 1.0) < 1e-12);

    c = case_c_min_time();
    assert(abs(c.p.A - 0.0) < 1e-12 && abs(c.p.B - 0.0) < 1e-12 && abs(c.p.C - 1.0) < 1e-12);
end
