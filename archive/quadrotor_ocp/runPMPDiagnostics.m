function [pmp, pmpSummary] = runPMPDiagnostics(problem, solution, lambda_struct, diagSpec, opts)
%==========================================================================
% runPMPDiagnostics.m  (DIAGNOSTICS MODULE)
%==========================================================================
% Responsibility:
%   Given a solved trajectory + fmincon multipliers, compute PMP/KKT signals:
%     - costate-like sequence from defect multipliers
%     - Hamiltonian trace
%     - dH/du, d2H/du2, d/dt(dH/du) via problem.pmp_derivatives()
%
% Inputs:
%   problem, solution (states, controls, time_vec)
%   lambda_struct from fmincon
%   diagSpec: explicit mapping assumptions (avoid hidden coupling)
%       diagSpec.defectMultiplierLength = ns*N
%       diagSpec.ceqLayout = "defects_then_terminal"
%   opts (optional):
%       opts.smooth_costates = true/false
%       opts.smooth_window_frac = 1/20
%
% IMPORTANT ASSUMPTION (explicit):
%   solver builds ceq as: [defects(:); terminal_eq(:)]
%   thus eqnonlin multipliers begin with defect multipliers.
%==========================================================================

    if nargin < 4 || isempty(diagSpec)
        error('diagSpec is required (pass defectMultiplierLength, ceqLayout).');
    end
    if nargin < 5, opts = struct(); end
    opts = with_defaults(opts, struct( ...
        'smooth_costates', true, ...
        'smooth_window_frac', 1/20 ...
    ));

    if ~isfield(problem,'pmp_derivatives') || isempty(problem.pmp_derivatives)
        error('problem.pmp_derivatives is empty; PMP diagnostics not available.');
    end

    X = solution.states;
    U = solution.controls;
    t = solution.time_vec(:)';
    T = solution.T_final;

    ns = problem.num_states;
    N  = problem.time.N;

    if diagSpec.defectMultiplierLength ~= ns*N
        warning('diagSpec.defectMultiplierLength (%d) != ns*N (%d). Proceeding anyway.', ...
            diagSpec.defectMultiplierLength, ns*N);
    end

    % --- Extract multipliers
    if ~isfield(lambda_struct,'eqnonlin') || isempty(lambda_struct.eqnonlin)
        error('lambda_struct.eqnonlin missing/empty.');
    end
    eqmult = lambda_struct.eqnonlin(:);

    Ldef = diagSpec.defectMultiplierLength;
    if numel(eqmult) < Ldef
        error('eqnonlin multipliers too short (%d). Need at least %d.', numel(eqmult), Ldef);
    end

    defect_mult = eqmult(1:Ldef);          % assumes defects first
    MU = reshape(defect_mult, ns, N);      % ns x N

    % --- Costate-like mapping (documented heuristic)
    Lambda_raw = -[MU, MU(:,end)];         % ns x (N+1)

    % --- Optional smoothing (readability)
    Lambda = Lambda_raw;
    if opts.smooth_costates
        win = max(3, round((N+1)*opts.smooth_window_frac));
        for i = 1:ns
            Lambda(i,:) = smoothdata(Lambda(i,:), 'gaussian', win);
        end
    end

    % --- Hamiltonian
    h = T/N;

    f_mat = problem.rhs_collocation(X, U, problem.params);  % ns x (N+1)

    Udot = zeros(size(U));
    for i = 1:problem.num_controls
        Udot(i,:) = gradient(U(i,:), h);
    end

    L_vec = problem.cost_integrand(X, U, Udot, problem.params, problem.cost_weights);
    L_vec = reshape(L_vec, 1, []);
    H_vec = L_vec + sum(Lambda_raw .* f_mat, 1);

    % --- Derivatives (stationarity checks)
    dH_du_all    = [];
    d2H_du2_all  = [];
    dH_du_dt_all = [];

    try
        [dH_du_all, d2H_du2_all, dH_du_dt_all] = problem.pmp_derivatives( ...
            X, U, Lambda_raw, problem.params, problem.cost_weights, t);
    catch ME_raw
        warning('pmp_derivatives failed with raw costates: %s. Retrying with smoothed costates.', ME_raw.message);
        [dH_du_all, d2H_du2_all, dH_du_dt_all] = problem.pmp_derivatives( ...
            X, U, Lambda, problem.params, problem.cost_weights, t);
    end

    % --- Package outputs
    pmp = struct();
    pmp.costates_raw = Lambda_raw;
    pmp.costates     = Lambda;
    pmp.hamiltonian  = H_vec;
    pmp.dH_du        = dH_du_all;
    pmp.d2H_du2      = d2H_du2_all;
    pmp.dH_du_dt     = dH_du_dt_all;

    pmpSummary = struct();
    pmpSummary.status = "ok";
    pmpSummary.H_mean = mean(H_vec);
    pmpSummary.H_std  = std(H_vec);

    if ~isempty(dH_du_all)
        pmpSummary.max_abs_dH_du = max(abs(dH_du_all(:)));
    else
        pmpSummary.max_abs_dH_du = NaN;
    end
end

function s = with_defaults(s, d)
    f = fieldnames(d);
    for k = 1:numel(f)
        if ~isfield(s, f{k}) || isempty(s.(f{k}))
            s.(f{k}) = d.(f{k});
        end
    end
end
