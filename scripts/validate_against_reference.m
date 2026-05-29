function report = validate_against_reference(summaries, varargin)
%VALIDATE_AGAINST_REFERENCE Compare computed metrics to paper Table 3.
%
% Usage:
%   report = validate_against_reference(summaries, 'verbose', true);
%
% If results/reference/canonical_outputs.mat exists, also compares trajectories
% at a coarse level (terminal tumor, T, exposures).

    p = inputParser;
    addParameter(p, 'verbose', true, @islogical);
    parse(p, varargin{:});
    verbose = p.Results.verbose;

    ref = reference_metrics();
    repro = reproducibility_settings();

    cases = {'case_a', 'case_b', 'case_c'};
    ref_cases = {ref.case_a, ref.case_b, ref.case_c};

    report = struct();
    report.all_passed = true;
    report.rows = struct([]);

    if verbose
        fprintf('\n--- Reference Validation (Paper Table 3) ---\n');
        fprintf('%-10s %-8s %-12s %-12s %-12s %-8s\n', ...
            'Case', 'Metric', 'Computed', 'Reference', 'Abs Err', 'Pass');
    end

    for k = 1:numel(cases)
        s = summaries.(cases{k});
        r = ref_cases{k};

        metrics = { ...
            struct('name', 'T', 'computed', s.T_final, 'reference', r.T, ...
                   'tol', repro.metric_abs_tol_T); ...
            struct('name', 'Du', 'computed', s.total_chemo, 'reference', r.Du, ...
                   'tol', repro.metric_abs_tol_exposure); ...
            struct('name', 'Dv', 'computed', s.total_immuno, 'reference', r.Dv, ...
                   'tol', repro.metric_abs_tol_exposure); ...
            struct('name', 'x(T)', 'computed', s.terminal_tumor, 'reference', 1.0, ...
                   'tol', repro.metric_abs_tol_xf)};

        for m = 1:numel(metrics)
            met = metrics(m);
            err = abs(met.computed - met.reference);
            passed = err <= met.tol;
            report.all_passed = report.all_passed && passed;

            row = struct();
            row.case = cases{k};
            row.metric = met.name;
            row.computed = met.computed;
            row.reference = met.reference;
            row.abs_error = err;
            row.tolerance = met.tol;
            row.passed = passed;
            report.rows(end + 1) = row; %#ok<AGROW>

            if verbose
                if passed
                    ps = 'PASS';
                else
                    ps = 'FAIL';
                end
                fprintf('%-10s %-8s %12.4f %12.4f %12.4f %8s\n', ...
                    cases{k}, met.name, met.computed, met.reference, err, ps);
            end
        end
    end

    repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    canon_file = fullfile(repo_root, 'results', 'reference', 'canonical_outputs.mat');
    report.has_frozen_reference = isfile(canon_file);
    if report.has_frozen_reference
        loaded = load(canon_file, 'reference');
        report.frozen_reference_date = loaded.reference.created;
        if verbose
            fprintf('\nFrozen reference found (%s). Metric comparison above is primary check.\n', ...
                loaded.reference.created);
        end
    elseif verbose
        fprintf('\nNo frozen canonical_outputs.mat yet. Run reproduce_paper(''freeze_reference'', true).\n');
    end

    if verbose
        if report.all_passed
            fprintf('Overall: PASSED\n');
        else
            fprintf('Overall: NEEDS REVIEW\n');
        end
    end
end
