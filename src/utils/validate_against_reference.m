function report = validate_against_reference(summaries, varargin)
%VALIDATE_AGAINST_REFERENCE Compare metrics to published Table 3 values.

    p = inputParser;
    addParameter(p, 'verbose', true, @islogical);
    parse(p, varargin{:});
    verbose = p.Results.verbose;

    ref = reference_metrics();
    repro = reproducibility_settings();
    cases = {'case_a', 'case_b', 'case_c'};
    ref_cases = {ref.case_a, ref.case_b, ref.case_c};

    report = struct('all_passed', true, 'rows', struct([]));

    if verbose
        fprintf('\n--- Table 3 check ---\n');
    end

    for k = 1:numel(cases)
        if ~isfield(summaries, cases{k})
            continue;
        end
        s = summaries.(cases{k});
        r = ref_cases{k};
        checks = {
            'T', s.T_final, r.T, repro.metric_abs_tol_T;
            'Du', s.total_chemo, r.Du, repro.metric_abs_tol_exposure;
            'Dv', s.total_immuno, r.Dv, repro.metric_abs_tol_exposure;
            'x(T)', s.terminal_tumor, 1.0, repro.metric_abs_tol_xf};

        for m = 1:size(checks, 1)
            err = abs(checks{m, 2} - checks{m, 3});
            passed = err <= checks{m, 4};
            report.all_passed = report.all_passed && passed;
            if verbose
                if passed, mark = 'ok'; else, mark = 'check'; end
                fprintf('  %s %-6s  computed=%.3f  ref=%.3f  %s\n', ...
                    cases{k}, checks{m, 1}, checks{m, 2}, checks{m, 3}, mark);
            end
        end
    end
end
