function test_smoke_overrides()
%TEST_SMOKE_OVERRIDES Smoke path shrinks grid/swarm and flags non-paper.

    cfg = apply_reproducibility_settings(case_a_balanced());
    smoke = apply_smoke_overrides(cfg);
    repro = reproducibility_settings();

    assert(smoke.N == repro.smoke_grid_N);
    assert(smoke.qpso.swarmSize == repro.smoke_qpso_swarm_size);
    assert(smoke.qpso.maxIter == repro.smoke_qpso_max_iter);
    assert(smoke.smoke == true);
    assert(contains(lower(smoke.case_name), 'smoke'));
    % Paper case weights must still be intact after smoke overrides
    assert(abs(smoke.p.A - 1.0) < 1e-12 && abs(smoke.p.B - 0.3) < 1e-12);
end
