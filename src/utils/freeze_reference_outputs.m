function freeze_reference_outputs(solutions, summaries, repo_root)
%FREEZE_REFERENCE_OUTPUTS Save canonical trajectories and metrics for validation.

    ref_dir = fullfile(repo_root, 'results', 'reference');
    if ~exist(ref_dir, 'dir')
        mkdir(ref_dir);
    end

    reference = struct();
    reference.created = datestr(now, 31);
    reference.repro = reproducibility_settings();
    reference.solutions = solutions;
    reference.summaries = summaries;
    reference.metrics = reference_metrics();

    save(fullfile(ref_dir, 'canonical_outputs.mat'), 'reference', '-v7.3');

    table3 = generate_table3(summaries, fullfile(ref_dir, 'table3_computed.csv'));

    meta = struct();
    meta.matlab_version = version;
    meta.reproducibility = reference.repro;
    meta.table3 = table3;
    save(fullfile(ref_dir, 'freeze_metadata.mat'), 'meta');

    fprintf('Frozen reference outputs under %s\n', ref_dir);
end
