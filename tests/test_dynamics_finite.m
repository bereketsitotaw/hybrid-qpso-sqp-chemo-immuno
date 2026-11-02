function test_dynamics_finite()
%TEST_DYNAMICS_FINITE One RHS evaluation stays finite under MM saturation.

    p = get_model_parameters();
    t = 0;
    w = [p.x0; p.y0; p.z0];
    time_vec = [0, 1];
    u_vec = [0.5, 0.5];
    v_vec = [0.25, 0.25];
    dw = dynamics(t, w, u_vec, v_vec, p, time_vec);
    assert(all(isfinite(dw)), 'dynamics produced non-finite derivatives');
    assert(numel(dw) == 3);
end
