function summaries = run_cases(case_id)
%RUN_CASES Solve and plot chemo-immunotherapy OCP case studies.

    repro = reproducibility_settings();
    rng(repro.rng_seed, repro.rng_generator);

    repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    results_dir = fullfile(repo_root, 'results');
    figures_dir = fullfile(results_dir, 'figures');
    if ~exist(figures_dir, 'dir')
        mkdir(figures_dir);
    end

    all_cases = struct( ...
        'a', @case_a_balanced, ...
        'b', @case_b_asymmetric, ...
        'c', @case_c_min_time);

    tags = fieldnames(all_cases);
    if strcmpi(case_id, 'all')
        run_list = tags;
    else
        case_id = lower(char(case_id));
        if ~isfield(all_cases, case_id)
            error('Unknown case "%s". Use ''a'', ''b'', ''c'', or omit for all.', case_id);
        end
        run_list = {case_id};
    end

    summaries = struct();
    fprintf('Hybrid QPSO-SQP | rng(%d) | MATLAB %s\n\n', repro.rng_seed, version);

    for k = 1:numel(run_list)
        tag = run_list{k};
        cfg = apply_reproducibility_settings(all_cases.(tag)());
        cfg.qpso.verbose = true;

        fprintf('--- %s ---\n', cfg.case_name);
        [solution, summary] = solve_chemo_immuno_ocp(cfg);
        print_performance_summary(solution, summary);

        summaries.(['case_' tag]) = summary;
        save_solution(solution, summary, fullfile(results_dir, ['case_' tag '_solution.mat']));
        plot_chemo_immuno_results(solution, summary, figures_dir);
    end

    if numel(run_list) == 3
        generate_table3(summaries, fullfile(results_dir, 'table3.csv'));
        validate_against_reference(summaries, 'verbose', true);
    end

    fprintf('\nDone. Results in: %s\n', results_dir);
end
