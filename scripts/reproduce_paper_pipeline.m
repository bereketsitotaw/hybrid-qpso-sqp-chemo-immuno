function report = reproduce_paper_pipeline(varargin)
%REPRODUCE_PAPER_PIPELINE One-command reproduction of Section 5 + Table 3.
%
%   reproduce_paper_pipeline()
%   reproduce_paper_pipeline('freeze_reference', true)

    p = inputParser;
    addParameter(p, 'freeze_reference', false, @islogical);
    addParameter(p, 'export_figures', true, @islogical);
    addParameter(p, 'run_sensitivity', false, @islogical);
    parse(p, varargin{:});
    opts = p.Results;

    repo_root = init_project();
    repro = reproducibility_settings();
    rng(repro.rng_seed, repro.rng_generator);

    fprintf('========================================\n');
    fprintf('PAPER REPRODUCTION PIPELINE\n');
    fprintf('MATLAB %s | rng(%d,''twister'')\n', version, repro.rng_seed);
    fprintf('========================================\n\n');

    case_fns = {@case_a_balanced, @case_b_asymmetric, @case_c_min_time};
    tags = {'case_a', 'case_b', 'case_c'};
    fig_tags = {'fig1_caseA', 'fig3_caseB', 'fig4_caseC'};

    results_dir = fullfile(repo_root, 'results');
    figures_dir = fullfile(results_dir, 'figures');
    tables_dir = fullfile(results_dir, 'tables');

    solutions = struct();
    summaries = struct();

    for k = 1:numel(case_fns)
        cfg = case_fns{k}();
        cfg = apply_reproducibility_settings(cfg);
        cfg.qpso.verbose = (k == 1);

        fprintf('\n>>> %s\n', cfg.case_name);
        [sol, summ] = solve_chemo_immuno_ocp(cfg);
        print_performance_summary(sol, summ);

        solutions.(tags{k}) = sol;
        summaries.(tags{k}) = summ;

        save_solution(sol, summ, fullfile(results_dir, [tags{k} '_solution.mat']));

        if opts.export_figures
            export_publication_figures(sol, fig_tags{k}, figures_dir, {'png', 'eps'});
        end
    end

    generate_table3(summaries, fullfile(tables_dir, 'table3.csv'));
    ref_dir = fullfile(results_dir, 'reference');
    if ~exist(ref_dir, 'dir')
        mkdir(ref_dir);
    end
    generate_table3(summaries, fullfile(ref_dir, 'table3_computed.csv'));

    report = validate_against_reference(summaries, 'verbose', true);

    if opts.freeze_reference
        freeze_reference_outputs(solutions, summaries, repo_root);
    end

    if opts.run_sensitivity
        fprintf('\n>>> Sensitivity analysis (Section 5.4, reduced grid)\n');
        run(fullfile(repo_root, 'scripts', 'run_sensitivity.m'));
    end

    fprintf('\n========================================\n');
    fprintf('REPRODUCTION COMPLETE\n');
    if report.all_passed
        fprintf('Reference validation: PASSED (within published tolerances)\n');
    else
        fprintf('Reference validation: REVIEW (see metrics above)\n');
    end
    fprintf('Outputs: %s\n', results_dir);
    fprintf('========================================\n');
end
