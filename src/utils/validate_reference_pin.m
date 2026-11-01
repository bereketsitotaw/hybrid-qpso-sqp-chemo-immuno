function report = validate_reference_pin(varargin)
%VALIDATE_REFERENCE_PIN Self-check that Table 3 pin matches published numbers.
%
% Does not require a solve or a .mat file — validates configs/reference_metrics.m
% against the hard-coded DCDS-B Table 3 values.

    p = inputParser;
    addParameter(p, 'verbose', true, @islogical);
    parse(p, varargin{:});
    verbose = p.Results.verbose;

    ref = reference_metrics();
    expected = struct( ...
        'case_a', struct('T', 4.917, 'Du', 2.771, 'Dv', 1.785), ...
        'case_b', struct('T', 4.723, 'Du', 2.528, 'Dv', 4.723), ...
        'case_c', struct('T', 4.417, 'Du', 3.401, 'Dv', 4.417));

    report = struct('all_passed', true);

    if verbose
        fprintf('\n--- Reference pin (Table 3) self-check ---\n');
    end

    cases = {'case_a', 'case_b', 'case_c'};
    for k = 1:numel(cases)
        c = cases{k};
        for f = {'T', 'Du', 'Dv'}
            got = ref.(c).(f{1});
            want = expected.(c).(f{1});
            ok = abs(got - want) <= 1e-9;
            report.all_passed = report.all_passed && ok;
            if verbose
                if ok, mark = 'ok'; else, mark = 'FAIL'; end
                fprintf('  %s %-2s  pin=%.3f  paper=%.3f  %s\n', c, f{1}, got, want, mark);
            end
        end
    end

    red_a = ref.chemo_reduction_vs_c.case_a_pct;
    red_b = ref.chemo_reduction_vs_c.case_b_pct;
    ok_a = abs(red_a - (-18.5)) <= 1e-9;
    ok_b = abs(red_b - (-25.7)) <= 1e-9;
    report.all_passed = report.all_passed && ok_a && ok_b;
    if verbose
        fprintf('  chemo reduction A=%.1f%% (paper -18.5)  B=%.1f%% (paper -25.7)\n', red_a, red_b);
        if report.all_passed
            fprintf('  pin OK\n');
        else
            fprintf('  pin FAILED\n');
        end
    end
end
