function summaries = run_cases(case_id)
%RUN_CASES Solve and plot chemo-immunotherapy OCP case studies.

    repro = reproducibility_settings();
    rng(repro.rng_seed, repro.rng_generator);

    repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    results_dir = fullfile(repo_root, 'results');
    generated_dir = fullfile(results_dir, 'generated');
    figures_dir = fullfile(results_dir, 'figures');
    if ~exist(figures_dir, 'dir')
        mkdir(figures_dir);
    end
    if ~exist(generated_dir, 'dir')
        mkdir(generated_dir);
    end

    smoke = false;
    if nargin < 1 || isempty(case_id)
        case_id = 'all';
    end
    case_id = lower(char(case_id));
    if strcmp(case_id, 'smoke')
        smoke = true;
        case_id = 'a';
    end

    all_cases = struct( ...
        'a', @case_a_balanced, ...
        'b', @case_b_asymmetric, ...
        'c', @case_c_min_time);

    tags = fieldnames(all_cases);
    if strcmpi(case_id, 'all')
        run_list = tags;
    else
        if ~isfield(all_cases, case_id)
            error('Unknown case "%s". Use ''a'', ''b'', ''c'', ''smoke'', or omit for all.', case_id);
        end
        run_list = {case_id};
    end

    summaries = struct();
    if smoke
        fprintf('Hybrid QPSO-SQP | SMOKE (non-paper) | rng(%d) | MATLAB %s\n\n', ...
            repro.rng_seed, version);
    else
        fprintf('Hybrid QPSO-SQP | publication defaults | rng(%d) | MATLAB %s\n\n', ...
            repro.rng_seed, version);
    end

    for k = 1:numel(run_list)
        tag = run_list{k};
        cfg = all_cases.(tag)();
        assert_publication_defaults(cfg);
        cfg = apply_reproducibility_settings(cfg);
        if smoke
            cfg = apply_smoke_overrides(cfg);
        else
            cfg.smoke = false;
            cfg.qpso.verbose = true;
        end

        fprintf('--- %s ---\n', cfg.case_name);
        fprintf('  (A,B,C)=(%.3g,%.3g,%.3g) | N=%d | swarm=%d | maxIter=%d\n', ...
            cfg.p.A, cfg.p.B, cfg.p.C, cfg.N, cfg.qpso.swarmSize, cfg.qpso.maxIter);

        [solution, summary] = solve_chemo_immuno_ocp(cfg);
        summary.smoke = smoke;
        summary.case_tag = tag;
        print_performance_summary(solution, summary);

        summaries.(['case_' tag]) = summary;
        mat_name = ['case_' tag '_solution.mat'];
        save_solution(solution, summary, fullfile(generated_dir, mat_name));
        % Legacy mirror for older scripts / README paths
        save_solution(solution, summary, fullfile(results_dir, mat_name));
        if ~smoke
            plot_chemo_immuno_results(solution, summary, figures_dir);
        end
    end

    if ~smoke && numel(run_list) == 3
        generate_table3(summaries, fullfile(results_dir, 'table3.csv'));
        validate_against_reference(summaries, 'verbose', true);
    end

    fprintf('\nDone. Generated: %s\n', generated_dir);
end
