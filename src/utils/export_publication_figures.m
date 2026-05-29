function export_publication_figures(solution, tag, out_dir, formats)
%EXPORT_PUBLICATION_FIGURES Save state/control and PMP panels with paper names.
%
% tag examples: 'fig1_caseA', 'fig2_pmp_caseA', 'fig3_caseB'

    if nargin < 4 || isempty(formats)
        formats = {'png', 'eps'};
    end

    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end

    plot_chemo_immuno_results(solution, struct('case_name', solution.case_name), '');

    figs = findall(0, 'Type', 'figure');
    figs = flipud(figs);

    if numel(figs) >= 2
        state_fig = figs(2);
        pmp_fig = figs(1);
    elseif numel(figs) == 1
        state_fig = figs(1);
        pmp_fig = [];
    else
        return;
    end

    base_state = fullfile(out_dir, [tag '_states_controls']);
    base_pmp = fullfile(out_dir, [tag '_pmp']);

    for f = 1:numel(formats)
        fmt = lower(formats{f});
        switch fmt
            case 'png'
                exportgraphics(state_fig, [base_state '.png'], 'Resolution', 300);
                if ~isempty(pmp_fig)
                    exportgraphics(pmp_fig, [base_pmp '.png'], 'Resolution', 300);
                end
            case 'eps'
                exportgraphics(state_fig, [base_state '.eps'], 'ContentType', 'vector');
                if ~isempty(pmp_fig)
                    exportgraphics(pmp_fig, [base_pmp '.eps'], 'ContentType', 'vector');
                end
            otherwise
                warning('Unsupported format: %s', fmt);
        end
    end

    close(state_fig);
    if ~isempty(pmp_fig)
        close(pmp_fig);
    end
end
