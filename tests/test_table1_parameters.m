function test_table1_parameters()
%TEST_TABLE1_PARAMETERS Table 1 biological parameters match DCDS-B.

    p = get_model_parameters();
    assert(abs(p.xi - 0.5618) < 1e-12);
    assert(abs(p.x_inf - 780) < 1e-12);
    assert(abs(p.a - 0.00726) < 1e-12);
    assert(abs(p.b - 0.00264) < 1e-12);
    assert(abs(p.gamma - 0.1181) < 1e-12);
    assert(abs(p.delta - 0.3745) < 1e-12);
    assert(abs(p.x0 - 650) < 1e-12);
    assert(abs(p.y0 - 0.1) < 1e-12);
    assert(abs(p.z0 - 400) < 1e-12);
    assert(abs(p.xf - 1.0) < 1e-12);
    assert(abs(p.umax - 1.0) < 1e-12);
    assert(abs(p.vmax - 1.0) < 1e-12);
end
