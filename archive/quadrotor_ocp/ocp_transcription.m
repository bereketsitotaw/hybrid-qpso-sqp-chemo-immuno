%==========================================================================
% ocp_transcription.m
%==========================================================================
% Returns a "transcription interface" struct with function handles:
%   tr.parse_z(z) -> [X,U,T]
%   tr.objective(z)
%   tr.constraints(z) -> [c,ceq]  (includes DENSE inequalities!)
%   tr.rollout(x0,U,T)
%   tr.metrics(X,U,T) -> standardized reporting metrics incl dense checks
%
% Dense inequality checks (between-node) are evaluated via Hermite sampling
% on [x,y] using [vx,vy] when available.
%==========================================================================
function tr = ocp_transcription(problem)

    tr.parse_z     = @(z) parse_z_impl(z, problem);
    tr.objective   = @(z) objective_impl(z, problem);
    tr.constraints = @(z) constraints_impl(z, problem);
    tr.rollout     = @(x0,U,T) rollout_trap_heun_impl(x0, U, T, problem);
    tr.metrics     = @(X,U,T) metrics_impl(X, U, T, problem);
end

%==========================================================================
% parse z
%==========================================================================
function [states, controls, T] = parse_z_impl(z, problem)
    N  = problem.time.N;
    ns = problem.num_states;
    nc = problem.num_controls;

    idx_states_end = ns*(N+1);
    states = reshape(z(1:idx_states_end), ns, N+1);

    idx_controls_end = idx_states_end + nc*(N+1);
    controls = reshape(z(idx_states_end+1:idx_controls_end), nc, N+1);

    if isfield(problem.time,'T_fixed') && ~isempty(problem.time.T_fixed)
        T = problem.time.T_fixed;
    else
        T = z(end);
    end
end

%==========================================================================
% objective
%==========================================================================
function J = objective_impl(z, problem)
    [X, U, T] = parse_z_impl(z, problem);
    N = problem.time.N;
    t = linspace(0, T, N+1);
    h = T/N;

    Udot = zeros(size(U));
    for i = 1:problem.num_controls
        Udot(i,:) = gradient(U(i,:), h);
    end

    L = problem.cost_integrand(X, U, Udot, problem.params, problem.cost_weights);
    L = reshape(L,1,[]);
    J = trapz(t, L);

    if isfield(problem,'cost_terminal') && ~isempty(problem.cost_terminal)
        J = J + problem.cost_terminal(X(:,end), problem.params, problem.cost_weights);
    end
end

%==========================================================================
% constraints: defects + terminal eq + (node ineq if provided) + dense ineq
%==========================================================================
function [c, ceq] = constraints_impl(z, problem)
    [X, U, T] = parse_z_impl(z, problem);

    N = problem.time.N;
    h = T/N;
    t = linspace(0, T, N+1);

    % Explicit trapezoid (Heun) defects (kept consistent with rollout)
    defects = zeros(problem.num_states, N);
    for k = 1:N
        xk  = X(:,k);
        uk  = U(:,k);
        uk1 = U(:,k+1);

        fk   = rhs_at_node_impl(xk, uk, problem);
        xpred = xk + h*fk;
        fk1p = rhs_at_node_impl(xpred, uk1, problem);

        defects(:,k) = X(:,k+1) - X(:,k) - (h/2)*(fk + fk1p);
    end
    ceq_def = defects(:);

    % terminal equalities (respect NaNs)
    xtgt = problem.bounds.terminal_states(:);
    xT   = X(:,end);
    mask = ~isnan(xtgt);
    ceq_T = xT(mask) - xtgt(mask);

    ceq = [ceq_def; ceq_T];

    % user/node inequalities (optional)
    c_node = [];
    if isfield(problem,'inequality_constraints') && ~isempty(problem.inequality_constraints)
        c_node = problem.inequality_constraints(X, U, t, problem.params, problem.bounds);
        c_node = c_node(:);
    end

    % dense inequalities (hard safety/implementability)
    [c_dense, ~] = eval_ineq_dense_impl(X, U, T, problem);

    c = [c_node; c_dense];
    c = c(:);
end

%==========================================================================
% rollout (Heun / explicit trapezoid) — consistent with defects
%==========================================================================
function X = rollout_trap_heun_impl(x0, U, T, problem)
    N  = problem.time.N;
    ns = problem.num_states;
    h  = T/N;

    X = zeros(ns, N+1);
    X(:,1) = x0;

    for k = 1:N
        xk  = X(:,k);
        uk  = U(:,k);
        uk1 = U(:,k+1);

        fk   = rhs_at_node_impl(xk, uk, problem);
        if any(~isfinite(fk)), X(:) = NaN; return; end

        xpred = xk + h*fk;

        fk1p = rhs_at_node_impl(xpred, uk1, problem);
        if any(~isfinite(fk1p)), X(:) = NaN; return; end

        X(:,k+1) = xk + (h/2)*(fk + fk1p);
        if any(~isfinite(X(:,k+1))), X(:) = NaN; return; end
    end
end

%==========================================================================
% RHS at a node (fast path if problem.rhs exists; fallback to rhs_collocation)
%==========================================================================
function fk = rhs_at_node_impl(x, u, problem)
    if isfield(problem,'rhs') && ~isempty(problem.rhs)
        fk = problem.rhs(x, u, problem.params);
        return;
    end

    % fallback: use rhs_collocation with 2 columns
    tmpX = [x, x];
    tmpU = [u, u];
    F = problem.rhs_collocation(tmpX, tmpU, problem.params);
    fk = F(:,1);
end

%==========================================================================
% metrics (paper-ready, consistent everywhere)
%==========================================================================
function mt = metrics_impl(X, U, T, problem)
    N = problem.time.N;
    t = linspace(0, T, N+1);
    h = T/N;

    % J
    Udot = zeros(size(U));
    for i = 1:problem.num_controls
        Udot(i,:) = gradient(U(i,:), h);
    end
    L = problem.cost_integrand(X, U, Udot, problem.params, problem.cost_weights);
    L = reshape(L,1,[]);
    J = trapz(t, L);
    if isfield(problem,'cost_terminal') && ~isempty(problem.cost_terminal)
        J = J + problem.cost_terminal(X(:,end), problem.params, problem.cost_weights);
    end

    % terminal error (NaN entries ignored)
    xtgt = problem.bounds.terminal_states(:);
    xT   = X(:,end);
    term_err2 = 0;
    for i = 1:problem.num_states
        if ~isnan(xtgt(i))
            e = xT(i) - xtgt(i);
            term_err2 = term_err2 + e^2;
        end
    end
    termErr = sqrt(term_err2);

    % dense inequalities
    [c_dense, info] = eval_ineq_dense_impl(X, U, T, problem);
    maxDense = max(0, max(c_dense, [], 'omitnan'));

    mt = struct();
    mt.J = J;
    mt.termErr = termErr;
    mt.maxDenseIneqViolation = maxDense;
    mt.minDenseClearance = info.minDenseClearance;
    mt.maxSpeedSlack = info.maxSpeedSlack;
    mt.maxRateSlack  = info.maxRateSlack;
end

%==========================================================================
% Dense inequality evaluator: obstacle/speed/map/rate (between-node safety)
% Convention: fmincon wants c <= 0. Positive => violation.
%
% Implemented for 2D quadrotor-like state layout by default:
%   X = [x;y;vx;vy; ...]
% Hermite interpolation on x,y using vx,vy. Speed from interpolated derivatives.
%==========================================================================
function [c_dense, info] = eval_ineq_dense_impl(X, U, T, problem)

    N  = problem.time.N;
    h  = T/N;

    % defaults / options
    M  = getfieldwithdefault(problem.options, 'dense_samples', 12);
    refine_factor = getfieldwithdefault(problem.options, 'dense_refine_factor', 2);
    refine_band   = getfieldwithdefault(problem.options, 'dense_refine_band', 0.25);

    % state layout
    pos_idx = [1 2];
    vel_idx = [3 4];
    if isfield(problem.options,'state_layout') && isstruct(problem.options.state_layout)
        if isfield(problem.options.state_layout,'pos_idx'), pos_idx = problem.options.state_layout.pos_idx; end
        if isfield(problem.options.state_layout,'vel_idx'), vel_idx = problem.options.state_layout.vel_idx; end
    end

    haveHermite = (problem.num_states >= max([pos_idx(:); vel_idx(:)]));

    % bounds (optional)
    v_max = [];
    if isfield(problem.bounds,'v_max') && ~isempty(problem.bounds.v_max)
        v_max = problem.bounds.v_max;
    end
    haveMap = isfield(problem.bounds,'xlim') && isfield(problem.bounds,'ylim') && ...
              numel(problem.bounds.xlim)==2 && numel(problem.bounds.ylim)==2;

    du_max = [];
    if isfield(problem.bounds,'du_max') && ~isempty(problem.bounds.du_max)
        du_max = problem.bounds.du_max(:);
        if numel(du_max)==1
            du_max = repmat(du_max, problem.num_controls, 1);
        end
    end

    % obstacle data (optional)
    obs = [];
    if isfield(problem.params,'obs') && ~isempty(problem.params.obs)
        obs = problem.params.obs;
    end

    safety_margin = 0;
    if isfield(problem.params,'safety_margin') && ~isempty(problem.params.safety_margin)
        safety_margin = problem.params.safety_margin;
    elseif isfield(problem.bounds,'safety_margin') && ~isempty(problem.bounds.safety_margin)
        safety_margin = problem.bounds.safety_margin;
    end

    % Collect dense constraints
    c_obs_all  = [];
    c_spd_all  = [];
    c_map_all  = [];

    minClear = inf;
    maxSpeedSlack = -inf;

    % time base
    t_nodes = linspace(0, T, N+1);

    for k = 1:N
        t0 = t_nodes(k);
        t1 = t_nodes(k+1);

        % base sampling in tau
        tau = linspace(0,1,M);
        t_s = t0 + tau*(t1-t0);

        % interpolate pos/vel on interval
        if haveHermite
            p0 = X(pos_idx, k);
            p1 = X(pos_idx, k+1);
            v0 = X(vel_idx, k);
            v1 = X(vel_idx, k+1);

            [p_s, v_s] = hermite_posvel(p0, v0, p1, v1, h, tau);
        else
            % fallback: linear interp on position only
            p0 = X(pos_idx, k);
            p1 = X(pos_idx, k+1);
            p_s = p0 + (p1-p0).*tau;
            v_s = nan(size(p_s));
        end

        % quick clearance estimate to decide refinement
        [c_obs, clear_min_local] = obstacle_ineq(p_s, t_s, obs, safety_margin);
        localMinClear = clear_min_local;
        minClear = min(minClear, localMinClear);

        margin = getfieldwithdefault(problem.options,'dense_margin',0.0);
        % refine when close to boundary (near-zero clearance band)
        if isfinite(localMinClear) && (localMinClear < (margin + refine_band*max(1e-6, abs(margin)+1)))
            tau2 = linspace(0,1,refine_factor*M);
            t_s  = t0 + tau2*(t1-t0);

            if haveHermite
                [p_s, v_s] = hermite_posvel(p0, v0, p1, v1, h, tau2);
            else
                p_s = p0 + (p1-p0).*tau2;
                v_s = nan(size(p_s));
            end
            [c_obs, clear_min_local] = obstacle_ineq(p_s, t_s, obs, safety_margin);
            minClear = min(minClear, clear_min_local);
        end

        c_obs_all = [c_obs_all; c_obs(:)]; %#ok<AGROW>

        % speed ineq (if available)
        if ~isempty(v_max) && all(isfinite(v_s(:)))
            spd = sqrt(sum(v_s.^2, 1));
            c_spd = spd - v_max;  % <= 0 is feasible
            c_spd_all = [c_spd_all; c_spd(:)]; %#ok<AGROW>
            maxSpeedSlack = max(maxSpeedSlack, max(spd - v_max)); %#ok<AGROW>
        end

        % map bounds (if provided)
        if haveMap
            xmin = problem.bounds.xlim(1); xmax = problem.bounds.xlim(2);
            ymin = problem.bounds.ylim(1); ymax = problem.bounds.ylim(2);

            x = p_s(1,:); y = p_s(2,:);
            c_map = [xmin - x; x - xmax; ymin - y; y - ymax]; % all <= 0 feasible
            c_map_all = [c_map_all; c_map(:)]; %#ok<AGROW>
        end
    end

    % control-rate constraints (node-based, but included in "dense set")
    c_rate_all = [];
    maxRateSlack = -inf;
    if ~isempty(du_max)
        dU = diff(U,1,2);               % nc x N
        c1 = dU - du_max;               % <=0
        c2 = -dU - du_max;              % <=0
        c_rate_all = [c1(:); c2(:)];
        maxRateSlack = max([max(c1(:)), max(c2(:))]);
    end

    % Combine all
    c_dense = [c_obs_all; c_spd_all; c_map_all; c_rate_all];
    c_dense = c_dense(:);

    % info for reporting
    info = struct();
    info.minDenseClearance = minClear;
    info.maxSpeedSlack = maxSpeedSlack;
    info.maxRateSlack  = maxRateSlack;
end

%==========================================================================
% Hermite interpolation for position+velocity
% p(tau) = h00 p0 + h10 h v0 + h01 p1 + h11 h v1
% v(tau) = dp/dt (time derivative)
%==========================================================================
function [p_s, v_s] = hermite_posvel(p0, v0, p1, v1, h, tau)

    tau = tau(:).';  % 1 x M

    h00 =  2*tau.^3 - 3*tau.^2 + 1;
    h10 =      tau.^3 - 2*tau.^2 + tau;
    h01 = -2*tau.^3 + 3*tau.^2;
    h11 =      tau.^3 - tau.^2;

    p_s = p0*h00 + (h*v0)*h10 + p1*h01 + (h*v1)*h11;

    % derivatives wrt tau
    dh00 = 6*tau.^2 - 6*tau;
    dh10 = 3*tau.^2 - 4*tau + 1;
    dh01 = -dh00;
    dh11 = 3*tau.^2 - 2*tau;

    dp_dtau = p0*dh00 + (h*v0)*dh10 + p1*dh01 + (h*v1)*dh11;
    v_s = (1/h) * dp_dtau;
end

%==========================================================================
% obstacle inequality: dist >= (R + safety_margin)  =>  (R+sm) - dist <= 0
% also returns min signed clearance (dist - (R+sm))
%==========================================================================
function [c_obs, min_clear] = obstacle_ineq(p_s, t_s, obs, safety_margin)

    if isempty(obs)
        c_obs = [];
        min_clear = inf;
        return;
    end

    x = p_s(1,:); y = p_s(2,:);
    min_clear = inf;

    c_list = [];

    for j = 1:numel(obs)
        pos0 = obs(j).pos0(:);
        if ~isfield(obs(j),'vel') || isempty(obs(j).vel)
            vel = [0;0];
        else
            vel = obs(j).vel(:);
        end
        R = obs(j).R;

        ox = pos0(1) + vel(1)*t_s;
        oy = pos0(2) + vel(2)*t_s;

        dist = sqrt((x-ox).^2 + (y-oy).^2);
        clear = dist - (R + safety_margin);
        min_clear = min(min_clear, min(clear));

        c_j = (R + safety_margin) - dist; % <=0 feasible
        c_list = [c_list; c_j(:)]; %#ok<AGROW>
    end

    c_obs = c_list;
end

function v = getfieldwithdefault(s, field, defaultVal)
    if isstruct(s) && isfield(s, field)
        v = s.(field);
    else
        v = defaultVal;
    end
end
