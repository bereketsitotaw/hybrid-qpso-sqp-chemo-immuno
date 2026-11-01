function test_publication_assert()
%TEST_PUBLICATION_ASSERT assert_publication_defaults accepts paper cases.

    cfg = case_a_balanced();
    assert_publication_defaults(cfg);
    assert_publication_defaults(case_b_asymmetric());
    assert_publication_defaults(case_c_min_time());
end
