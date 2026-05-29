%==========================================================================
% solveOCP.m  (UNIFIED SOLVER PRODUCT — FIXED + REFACTORED)
%==========================================================================
% Hybrid method:
%   (1) QPSO coarse search over LOW-DIMENSIONAL control knots
%   (2) fmincon (SQP) direct transcription refinement on full grid
%   (3) OPTIONAL PMP diagnostics via external module runPMPDiagnostics.m
%
% Key fixes in this version:
%   - Dense (between-node) inequality checking via Hermite sampling
%   - Same dense-check used in QPSO, SQP constraints, and post-check
%   - Feasibility-first QPSO scoring (maxViolDense prioritized)
%   - Optional hard control-rate bound |u_{k+1}-u_k| <= du_max
%   - Standardized success metrics + eval budget accounting
%
% Dependencies (new files):
%   - qpso_knots_search.m
%   - ocp_transcription.m
%==========================================================================

function [solution, summary] = solveOCP(problem)

    %% 0) Header
    fprintf('========================================\n');
    fprintf('Starting OCP Solver for: %s\n', problem.name);
    fprintf('========================================\n');

    %% 1) Validate problem struct
    fprintf('1. Validating problem and setting up...\n');

    required_fields = {'num_states','num_controls','time','bounds', ...
                       'rhs_collocation','cost_integrand','params'};
    for f = required_fields
        if ~isfield(problem, f{1})
            error('Problem definition missing required field: %s', f{1});
        end
    end

    if ~isfield(problem,'cost_weights') || isempty(problem.cost_weights)
        problem.cost_weights = struct();
    end
    if ~isfield(problem,'options') || isempty(problem.options)
        problem.options = struct();
    end

    % Solver options
    problem.options = with_defaults(problem.options, struct( ...
        'enablePMP', true, ...
        'feas_tol_for_pmp', 1e-3, ...
        'dense_samples', 12, ...              % base dense samples per interval
        'dense_refine_factor', 2, ...         % refinement multiplier
        'dense_refine_band', 0.25, ...        % refine when clearance within band*margin
        'success_tol_ineq', 1e-6, ...
        'success_tol_term', 1e-3 ...
    ));

    % Dims
    N  = problem.time.N;
    ns = problem.num_states;
    nc = problem.num_controls;

    assert(numel(problem.bounds.initial_states)==ns, 'initial_states must be ns x 1');
    assert(numel(problem.bounds.control_min)==nc,    'control_min must be nc x 1');
    assert(numel(problem.bounds.control_max)==nc,    'control_max must be nc x 1');
    if isfield(problem.bounds,'terminal_states')
        assert(numel(problem.bounds.terminal_states)==ns, 'terminal_states must be ns x 1');
    else
        error('Problem definition missing required field: bounds.terminal_states');
    end

    % Time handling
    if isfield(problem.time,'T_fixed') && ~isempty(problem.time.T_fixed)
        T_guess       = problem.time.T_fixed;
        is_time_fixed = true;
    else
        if ~isfield(problem.time,'T_guess') || isempty(problem.time.T_guess)
            error('Must provide problem.time.T_fixed or problem.time.T_guess');
        end
        T_guess       = problem.time.T_guess;
        is_time_fixed = false;
    end

    % Build transcription interface (objective/constraints/rollout/dense checks)
    tr = ocp_transcription(problem);

    %% 2) QPSO defaults
    if ~isfield(problem,'qpso') || isempty(problem.qpso)
        problem.qpso = struct();
    end
    qpso = with_defaults(problem.qpso, struct( ...
        'swarmSize',   120, ...
        'maxIter',     120, ...
        'nKnots',      10, ...
        'nStages',     3, ...
        'useParfor',   true, ...
        'smooth_guess',true, ...
        'beta0',       1.0, ...   % QPSO contraction annealing
        'betaF',       0.45, ...
        'noise0',      0.06, ...
        'noiseF',      0.01, ...
        'stallIters',  30, ...
        'restartFrac', 0.25, ...
        'adaptKnots',  true, ...
        'adaptKnotsStep', 6, ...
        'adaptClearEps',  1e-2 ...
    ));

    %% 3) QPSO coarse search (knots)
    fprintf('2. Starting QPSO coarse search (knots=%d, T=%.2f)...\n', qpso.nKnots, T_guess);
    qpso_tic = tic;

    [knotsBest, U_guess, qpsoStats] = qpso_knots_search(problem, T_guess, tr, qpso);

    qpso_time = toc(qpso_tic);
    fprintf('   -> QPSO phase complete in %.2f seconds.\n', qpso_time);

    %% 3b) optional smoothing of the node-wise guess
    if qpso.smooth_guess
        umin = problem.bounds.control_min(:);
        umax = problem.bounds.control_max(:);
        for i = 1:nc
            U_guess(i,:) = smoothdata(U_guess(i,:), 'gaussian', max(3, round((N+1)/12)));
            U_guess(i,:) = min(max(U_guess(i,:), umin(i)), umax(i));
        end
    end

    %% 4) fmincon refinement
    fprintf('3. Starting fmincon local refinement...\n');
    fmin_tic = tic;

    X_guess = tr.rollout(problem.bounds.initial_states(:), U_guess, T_guess);

    z0_states   = reshape(X_guess, 1, ns*(N+1));
    z0_controls = reshape(U_guess, 1, nc*(N+1));

    lb_states = -Inf(1, ns*(N+1));
    ub_states =  Inf(1, ns*(N+1));
    lb_states(1:ns) = problem.bounds.initial_states(:);
    ub_states(1:ns) = problem.bounds.initial_states(:);

    lb_controls_flat = reshape(repmat(problem.bounds.control_min(:)', N+1, 1)', 1, nc*(N+1));
    ub_controls_flat = reshape(repmat(problem.bounds.control_max(:)', N+1, 1)', 1, nc*(N+1));

    if is_time_fixed
        lb = [lb_states, lb_controls_flat]';
        ub = [ub_states, ub_controls_flat]';
        z0 = [z0_states, z0_controls]';
    else
        if isfield(problem.time,'T_bounds') && numel(problem.time.T_bounds)==2
            lb_T = problem.time.T_bounds(1);
            ub_T = problem.time.T_bounds(2);
        else
            lb_T = 0.1; ub_T = 50.0;
        end
        % clamp initial guess into bounds
        T_guess = min(max(T_guess, lb_T), ub_T);

        lb = [lb_states, lb_controls_flat, lb_T]';
        ub = [ub_states, ub_controls_flat, ub_T]';
        z0 = [z0_states, z0_controls, T_guess]';
    end

    options = optimoptions('fmincon', ...
        'Algorithm','sqp', ...
        'Display','iter', ...
        'MaxIterations', 5000, ...
        'MaxFunctionEvaluations', 3e6, ...
        'ConstraintTolerance', 1e-6, ...
        'OptimalityTolerance', 1e-6 ...
    );

    [z_opt, J_opt, exitflag, output, lambda_struct] = fmincon( ...
        @(z) tr.objective(z), ...
        z0, [], [], [], [], lb, ub, ...
        @(z) tr.constraints(z), ...
        options );

    fmincon_time = toc(fmin_tic);
    fprintf('   -> fmincon refinement complete in %.2f seconds. Final cost J = %.6f\n', ...
        fmincon_time, J_opt);

    %% 5) Package solution & summary
    fprintf('4. Packaging final solution...\n');

    [X_opt, U_opt, T_opt] = tr.parse_z(z_opt);
    t_vec = linspace(0, T_opt, N+1);

    solution = struct();
    solution.problem_name = problem.name;
    solution.time_vec     = t_vec;
    solution.states       = X_opt;
    solution.controls     = U_opt;
    solution.T_final      = T_opt;
    solution.cost         = J_opt;
    solution.pmp          = [];
    solution.pmp_status   = "not_run";

    % Compute standardized dense feasibility metrics
    metrics_opt = tr.metrics(X_opt, U_opt, T_opt);

    summary = struct();
    summary.problem_name   = problem.name;

    summary.qpso_time      = qpso_time;
    summary.fmincon_time   = fmincon_time;
    summary.total_time     = qpso_time + fmincon_time;

    summary.exitflag       = exitflag;
    summary.fmincon_output = output;
    summary.lambda_struct  = lambda_struct;

    % Budget accounting
    summary.qpso_evals     = qpsoStats.qpso_evals;
    summary.fmincon_evals  = getfieldwithdefault(output, 'funcCount', NaN);

    % QPSO best metrics (dense)
    summary.qpso_best_cost             = qpsoStats.bestCost;
    summary.qpso_best_maxDenseIneqViol = qpsoStats.bestMetrics.maxDenseIneqViolation;
    summary.qpso_best_minDenseClear    = qpsoStats.bestMetrics.minDenseClearance;
    summary.qpso_best_termErr          = qpsoStats.bestMetrics.termErr;

    % SQP feasibility (both node and dense)
    [c_node, ceq_node] = tr.constraints(z_opt);
    summary.max_node_ineq_violation = max(0, max(c_node(:), [], 'omitnan'));
    summary.max_node_eq_violation   = max(abs(ceq_node(:)), [], 'omitnan');
    summary.max_dense_ineq_violation = metrics_opt.maxDenseIneqViolation;
    summary.min_dense_clearance      = metrics_opt.minDenseClearance;
    summary.term_err                 = metrics_opt.termErr;
    summary.max_speed_slack          = metrics_opt.maxSpeedSlack;
    summary.max_rate_slack           = metrics_opt.maxRateSlack;

    % Success criterion (paper-ready)
    summary.success = (summary.max_dense_ineq_violation <= problem.options.success_tol_ineq) && ...
                      (summary.term_err <= problem.options.success_tol_term) && ...
                      (exitflag > 0);

    if isfield(output,'constrviolation')
        summary.constrviolation = output.constrviolation;
    else
        summary.constrviolation = max([summary.max_node_ineq_violation, summary.max_node_eq_violation]);
    end

    fprintf('   -> Node-feas: max|ceq|=%.3e | max(max(0,c))=%.3e\n', ...
        summary.max_node_eq_violation, summary.max_node_ineq_violation);
    fprintf('   -> Dense-feas: maxDense=%.3e | minClear=%.3e | termErr=%.3e | success=%d\n', ...
        summary.max_dense_ineq_violation, summary.min_dense_clearance, summary.term_err, summary.success);

    %% 6) Optional PMP hook (compute only)
    doPMP = problem.options.enablePMP ...
        && isfield(problem,'pmp_derivatives') && ~isempty(problem.pmp_derivatives) ...
        && exitflag > 0 ...
        && summary.max_dense_ineq_violation < problem.options.feas_tol_for_pmp ...
        && summary.max_node_eq_violation < problem.options.feas_tol_for_pmp;

    if doPMP
        fprintf('5. Running PMP diagnostics module...\n');

        diagSpec = struct();
        diagSpec.defectMultiplierLength = ns*N;          % critical slice length
        diagSpec.ceqLayout = "defects_then_terminal";    % assumption documentation

        try
            [pmp, pmpSummary] = runPMPDiagnostics(problem, solution, lambda_struct, diagSpec, struct());
            solution.pmp = pmp;
            solution.pmp_status = "ok";
            summary.pmpSummary = pmpSummary;
        catch ME
            warning('PMP module failed: %s', ME.message);
            solution.pmp_status = "failed";
            summary.pmpSummary = struct('status',"failed",'message',ME.message);
        end
    else
        if ~problem.options.enablePMP
            solution.pmp_status = "disabled";
        elseif exitflag <= 0
            solution.pmp_status = "skipped_nonconverged";
        else
            solution.pmp_status = "skipped_infeasible";
        end
    end

    fprintf('========================================\n');
    fprintf('Solver finished.\n');
    fprintf('========================================\n');
end

%==========================================================================
% tiny local helper (kept local intentionally)
%==========================================================================
function s = with_defaults(s, d)
    f = fieldnames(d);
    for k = 1:numel(f)
        if ~isfield(s, f{k}) || isempty(s.(f{k}))
            s.(f{k}) = d.(f{k});
        end
    end
end

function v = getfieldwithdefault(s, field, defaultVal)
    if isstruct(s) && isfield(s, field)
        v = s.(field);
    else
        v = defaultVal;
    end
end
