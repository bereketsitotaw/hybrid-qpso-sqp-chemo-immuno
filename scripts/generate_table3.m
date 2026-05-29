%GENERATE_TABLE3 Recompute and export Table 3 from saved or fresh solutions.
addpath(fileparts(mfilename('fullpath')));
repo_root = init_project();

mat_a = fullfile(repo_root, 'results', 'case_a_solution.mat');
mat_b = fullfile(repo_root, 'results', 'case_b_solution.mat');
mat_c = fullfile(repo_root, 'results', 'case_c_solution.mat');

if isfile(mat_a) && isfile(mat_b) && isfile(mat_c)
    sa = load(mat_a, 'summary');
    sb = load(mat_b, 'summary');
    sc = load(mat_c, 'summary');
    summaries = struct('case_a', sa.summary, 'case_b', sb.summary, 'case_c', sc.summary);
else
    fprintf('Saved solutions not found — running full reproduction first...\n');
    reproduce_paper_pipeline('export_figures', false);
    sa = load(mat_a, 'summary');
    sb = load(mat_b, 'summary');
    sc = load(mat_c, 'summary');
    summaries = struct('case_a', sa.summary, 'case_b', sb.summary, 'case_c', sc.summary);
end

tables_dir = fullfile(repo_root, 'results', 'tables');
generate_table3(summaries, fullfile(tables_dir, 'table3.csv'));
validate_against_reference(summaries, 'verbose', true);
