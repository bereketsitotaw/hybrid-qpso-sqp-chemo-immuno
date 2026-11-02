function results = run_smoke_tests()
%RUN_SMOKE_TESTS Lightweight MATLAB checks (no full paper solve required).
%
%   From repo root:
%     main('test')
%     run('tests/run_smoke_tests.m')

    here = fileparts(mfilename('fullpath'));
    root = fileparts(here);
    addpath(fullfile(root, 'src', 'chemo_immuno'));
    addpath(fullfile(root, 'src', 'utils'));
    addpath(fullfile(root, 'configs'));
    addpath(here);

    tests = { ...
        @test_table1_parameters, ...
        @test_case_weights, ...
        @test_reference_metrics_pin, ...
        @test_dynamics_finite, ...
        @test_publication_assert, ...
        @test_smoke_overrides ...
        };

    names = { ...
        'table1_parameters', ...
        'case_weights', ...
        'reference_metrics_pin', ...
        'dynamics_finite', ...
        'publication_assert', ...
        'smoke_overrides' ...
        };

    n = numel(tests);
    passed = false(n, 1);
    messages = strings(n, 1);

    fprintf('=== hybrid-qpso-sqp-chemo-immuno smoke tests ===\n');
    for i = 1:n
        try
            tests{i}();
            passed(i) = true;
            messages(i) = "ok";
            fprintf('  PASS  %s\n', names{i});
        catch ME
            passed(i) = false;
            messages(i) = string(ME.message);
            fprintf('  FAIL  %s — %s\n', names{i}, ME.message);
        end
    end

    results = struct('names', {names}, 'passed', passed, 'messages', messages);
    fprintf('=== %d/%d passed ===\n', nnz(passed), n);
    assert(all(passed), '%d smoke test(s) failed', nnz(~passed));
end
