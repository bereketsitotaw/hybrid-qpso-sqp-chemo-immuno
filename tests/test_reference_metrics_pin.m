function test_reference_metrics_pin()
%TEST_REFERENCE_METRICS_PIN configs/reference_metrics.m matches Table 3.

    report = validate_reference_pin('verbose', false);
    assert(report.all_passed, 'Table 3 pin mismatch in reference_metrics.m');
end
