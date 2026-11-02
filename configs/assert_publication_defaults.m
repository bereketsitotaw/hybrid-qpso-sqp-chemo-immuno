function assert_publication_defaults(cfg)
%ASSERT_PUBLICATION_DEFAULTS Lock Table 1 / case weights to DCDS-B values.
%
% Fails fast if model parameters or case (A,B,C,l) drift from
% doi:10.3934/dcdsb.2026072. Does not invent alternate weight profiles.

    if nargin < 1 || isempty(cfg)
        error('assert_publication_defaults:MissingCfg', 'cfg struct required.');
    end

    p = get_model_parameters();
    expected = struct( ...
        'xi', 0.5618, 'x_inf', 780, 'theta', 1.0, 'alpha', 1.0, ...
        'a', 0.00726, 'b', 0.00264, 'gamma', 0.1181, 'delta', 0.3745, ...
        'kappa', 2.0, 'nu', 1.0, 'sigma', 0.1, 'psi', 2.0, 'mu', 0.45, ...
        'x0', 650, 'y0', 0.1, 'z0', 400, 'xf', 1.0, ...
        'umax', 1.0, 'vmax', 1.0);

    fields = fieldnames(expected);
    for k = 1:numel(fields)
        f = fields{k};
        if abs(p.(f) - expected.(f)) > 1e-12
            error('assert_publication_defaults:Table1Drift', ...
                'Table 1 parameter %s drifted: got %g, expected %g.', ...
                f, p.(f), expected.(f));
        end
    end

    if ~isfield(cfg, 'p')
        error('assert_publication_defaults:MissingP', 'cfg.p required.');
    end

    % Case-level weights must match the case factory used to build cfg.
    if isfield(cfg, 'case_name')
        name = lower(cfg.case_name);
        if contains(name, 'balanced')
            expect_ABC(cfg.p, 1.0, 0.3, 1.0);
        elseif contains(name, 'asymmetric')
            expect_ABC(cfg.p, 1.0, 0.0, 1.0);
        elseif contains(name, 'minimum-time') || contains(name, 'min-time')
            expect_ABC(cfg.p, 0.0, 0.0, 1.0);
        end
    end

    for f = {'l1', 'l2', 'l3'}
        if ~isfield(cfg.p, f{1}) || abs(cfg.p.(f{1})) > 1e-12
            error('assert_publication_defaults:LagrangianDrift', ...
                'Paper Section 5 cases use ell=(0,0,0); %s invalid.', f{1});
        end
    end
end

function expect_ABC(p, A, B, C)
    if abs(p.A - A) > 1e-12 || abs(p.B - B) > 1e-12 || abs(p.C - C) > 1e-12
        error('assert_publication_defaults:CaseWeightDrift', ...
            'Case (A,B,C) drifted: got (%.4g,%.4g,%.4g), expected (%.4g,%.4g,%.4g).', ...
            p.A, p.B, p.C, A, B, C);
    end
end
