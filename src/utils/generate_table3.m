function table3 = generate_table3(summaries, out_csv)
%GENERATE_TABLE3 Build Table 3-style comparison and optionally export CSV.
%
% summaries: struct with fields case_a, case_b, case_c (each a summary struct)

    ref = reference_metrics();

    cases = {'case_a', 'case_b', 'case_c'};
    labels = {ref.case_a.name, ref.case_b.name, ref.case_c.name};

    T_vals = zeros(1, 3);
    Du_vals = zeros(1, 3);
    Dv_vals = zeros(1, 3);
    Du_paper = [ref.case_a.Du, ref.case_b.Du, ref.case_c.Du];

    for k = 1:3
        s = summaries.(cases{k});
        T_vals(k) = s.T_final;
        Du_vals(k) = s.total_chemo;
        Dv_vals(k) = s.total_immuno;
    end

    Du_c = ref.case_c.Du;
    chemo_reduction = 100 * (Du_vals - Du_c) / Du_c;

    Metric = [ ...
        "Time (T)"; ...
        "Chemo (Du)"; ...
        "Immuno (Dv)"; ...
        "Chemo reduction vs Case C (%)"];
    CaseA = [T_vals(1); Du_vals(1); Dv_vals(1); chemo_reduction(1)];
    CaseB = [T_vals(2); Du_vals(2); Dv_vals(2); chemo_reduction(2)];
    CaseC = [T_vals(3); Du_vals(3); Dv_vals(3); 0];

    table3 = table(Metric, CaseA, CaseB, CaseC);

    fprintf('\n--- Table 3: Comparative Performance Metrics ---\n');
    disp(table3);

    if nargin >= 2 && ~isempty(out_csv)
        out_dir = fileparts(out_csv);
        if ~exist(out_dir, 'dir')
            mkdir(out_dir);
        end
        writetable(table3, out_csv);
        fprintf('Saved Table 3 to %s\n', out_csv);
    end
end
