function test_dynamics()
%TEST_DYNAMICS Sanity check for Michaelis-Menten RHS at nominal ICs.

    p = get_model_parameters();
    p.l1 = 0; p.l2 = 0; p.l3 = 0; p.A = 1; p.B = 0.3; p.C = 1;

    u_vec = zeros(1, 101);
    v_vec = zeros(1, 101);
    time = linspace(0, 4, 101);
    w0 = [p.x0; p.y0; p.z0];

    dw = dynamics(0, w0, u_vec, v_vec, p, time);

    assert(numel(dw) == 3, 'Dynamics must return 3 derivatives.');
    assert(all(isfinite(dw)), 'Dynamics must be finite at t=0.');
    assert(dw(1) > 0, 'Tumor should grow initially without therapy.');

    fprintf('test_dynamics: PASSED\n');
end
