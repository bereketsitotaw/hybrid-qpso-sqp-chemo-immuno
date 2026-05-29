function problem = define_quadrotor_nav_problem(opts)
%==========================================================================
% define_quadrotor_nav_problem.m  (DROP-IN FIXED)
%==========================================================================
% 2D quadrotor nav in clutter:
% States:   [x; y; vx; vy]
% Controls: [ax; ay]
%
% Constraints (intended):
%   - speed: ||v|| <= v_max
%   - obstacles: dist^2 >= (R + safety_margin)^2
%   - (optional) control rate bound: |u_{k+1}-u_k| <= du_max
%   - (optional) map bounds: xlim/ylim
%
% IMPORTANT (FIXED):
%   - Inequality constraints are evaluated on an UPSAMPLED grid inside each
%     interval using HERMITE interpolation (x,y from vx,vy) to reduce tunneling.
%   - Also populates problem.bounds.* fields (v_max, xlim, ylim, safety_margin)
%     so the refactored solver dense-check can enforce the same rules.
%==========================================================================

    if nargin < 1 || isempty(opts), opts = struct(); end
    fprintf('Loading quadrotor navigation problem...\n');

    opts.seed = getfield_default(opts,'seed',1);
    rng(opts.seed);

    % Dimensions
    problem.name         = 'Quadrotor Min-Time Nav (Clutter + Certified)';
    problem.num_states   = 4;
    problem.num_controls = 2;

    % Discretization
    problem.time.N = getfield_default(opts,'N',40);

    % Time settings
    problem.time.T_fixed  = [];
    problem.time.T_guess  = getfield_default(opts,'T_guess',6.0);
    problem.time.T_bounds = getfield_default(opts,'T_bounds',[0.2, 50.0]);

    % Bounds
    b = struct();
    b.initial_states  = [0; 0; 0; 0];
    b.terminal_states = [10; 10; 0; 0];

    a_max = getfield_default(opts,'a_max',2.0);
    b.control_min = [-a_max; -a_max];
    b.control_max = [ a_max;  a_max];

    % Optional control-rate bound at nodes: |u_{k+1}-u_k| <= du_max
    if isfield(opts,'du_max') && ~isempty(opts.du_max)
        if isscalar(opts.du_max), b.du_max = [opts.du_max; opts.du_max];
        else, b.du_max = opts.du_max(:);
        end
        assert(numel(b.du_max)==2,'du_max must be scalar or 2x1');
    else
        b.du_max = [];
    end

    % Parameters
    p = struct();
    p.v_max = getfield_default(opts,'max_velocity',5.0);
    p.safety_margin = getfield_default(opts,'safety_margin',0.40);

    p.map.xlim = getfield_default(opts,'xlim',[0,10]);
    p.map.ylim = getfield_default(opts,'ylim',[0,10]);

    num_obs    = getfield_default(opts,'num_obs',8);
    moving_obs = getfield_default(opts,'moving_obs',false);
    p.obs = make_clutter_obstacles(num_obs, moving_obs, p.map);

    % Legacy single obstacle compatibility (optional)
    p.obstacle = struct();
    p.obstacle.initial_pos = p.obs(1).pos0;
    p.obstacle.velocity    = p.obs(1).vel;
    p.obstacle.radius      = p.obs(1).R;

    % Upsample count inside each interval (for problem.inequality_constraints)
    p.constraint_upsample = max(1, round(getfield_default(opts,'constraint_upsample', 3)));

    problem.params = p;

    % ---- Compatibility: populate bounds fields used by solver dense-check ----
    % Speed / map bounds / safety margin duplicated in bounds for solver-side checks
    b.v_max         = p.v_max;
    b.safety_margin = p.safety_margin;
    b.xlim          = p.map.xlim;
    b.ylim          = p.map.ylim;

    problem.bounds = b;

    % Cost weights / objective mode
    w = struct();
    w.objective_mode = getfield_default(opts,'objective_mode','time_energy_L2');
    w.C     = getfield_default(opts,'C',1.0);
    w.wu    = getfield_default(opts,'wu',[1.0; 1.0]);
    w.wudot = getfield_default(opts,'wudot',[0.5; 0.5]);
    w.wu_reg_in_jerk = getfield_default(opts,'wu_reg_in_jerk',[1e-3; 1e-3]);

    % Mode -> time fixed/free
    switch w.objective_mode
        case 'time_energy_L2'
            if isfield(opts,'T_fixed') && ~isempty(opts.T_fixed)
                warning('time_energy_L2 uses FREE final time; ignoring opts.T_fixed.');
            end
            if w.C <= 0, error('Need C>0 for free-time well-posedness.'); end
            problem.time.T_fixed = [];
        case 'min_energy_L2'
            if ~isfield(opts,'T_fixed') || isempty(opts.T_fixed)
                error('min_energy_L2 requires opts.T_fixed.');
            end
            problem.time.T_fixed = opts.T_fixed;
        case 'min_jerk_L2'
            if ~isfield(opts,'T_fixed') || isempty(opts.T_fixed)
                error('min_jerk_L2 requires opts.T_fixed.');
            end
            problem.time.T_fixed = opts.T_fixed;
        otherwise
            error('Unknown objective_mode: %s', w.objective_mode);
    end
    problem.cost_weights = w;

    % QPSO params (solver will fill defaults too)
    qpso = struct();
    qpso.swarmSize = getfield_default(opts,'swarmSize',120);
    qpso.maxIter   = getfield_default(opts,'maxIter',120);
    qpso.nKnots    = getfield_default(opts,'nKnots',10);
    qpso.nStages   = getfield_default(opts,'nStages',3);
    problem.qpso   = qpso;

    % Options (PMP on/off + state layout for dense Hermite sampling)
    problem.options = struct();
    problem.options.enablePMP        = getfield_default(opts,'enablePMP',true);
    problem.options.feas_tol_for_pmp = getfield_default(opts,'feas_tol_for_pmp',1e-3);

    % State layout (used by solver-side dense checker if present)
    problem.options.state_layout = struct('pos_idx',[1 2], 'vel_idx',[3 4]);

    % If you want solver-side dense checker to run more/less samples, set here:
    if isfield(opts,'dense_samples') && ~isempty(opts.dense_samples)
        problem.options.dense_samples = opts.dense_samples;
    end
    if isfield(opts,'dense_refine_factor') && ~isempty(opts.dense_refine_factor)
        problem.options.dense_refine_factor = opts.dense_refine_factor;
    end
    if isfield(opts,'dense_refine_band') && ~isempty(opts.dense_refine_band)
        problem.options.dense_refine_band = opts.dense_refine_band;
    end

    % Function handles (problem only)
    problem.rhs_collocation        = @quadrotor_rhs_collocation;
    problem.cost_integrand         = @quadrotor_cost_integrand;
    problem.inequality_constraints = @quadrotor_ineq_constraints; % Hermite upsampled
    problem.generate_summary       = @quadrotor_summary;

    % PMP derivatives only for modes where L depends on u (not jerk unless augmented)
    if any(strcmp(w.objective_mode, {'time_energy_L2','min_energy_L2'}))
        problem.pmp_derivatives = @quadrotor_pmp_derivatives;
    else
        problem.pmp_derivatives = [];
    end
end

%==========================================================================
% Local functions (problem-specific only)
%==========================================================================

function obs = make_clutter_obstacles(num_obs, moving_obs, map)
    obs = repmat(struct('pos0',[0;0],'vel',[0;0],'R',1.0), 1, num_obs);

    Rmin = 0.6; Rmax = 1.2;
    start = [map.xlim(1); map.ylim(1)];
    goal  = [map.xlim(2); map.ylim(2)];

    for k = 1:num_obs
        placed = false;
        for tries = 1:400
            pos0 = [ map.xlim(1) + (map.xlim(2)-map.xlim(1))*rand;
                     map.ylim(1) + (map.ylim(2)-map.ylim(1))*rand ];
            R = Rmin + (Rmax-Rmin)*rand;

            if norm(pos0-start) < 2.0 || norm(pos0-goal) < 2.0
                continue;
            end

            ok = true;
            for j = 1:k-1
                if norm(pos0 - obs(j).pos0) < (R + obs(j).R + 0.4)
                    ok = false; break;
                end
            end
            if ~ok, continue; end

            obs(k).pos0 = pos0;
            obs(k).R    = R;

            if moving_obs
                obs(k).vel = 0.30 * (rand(2,1)-0.5);
            else
                obs(k).vel = [0;0];
            end

            placed = true;
            break;
        end

        if ~placed
            obs(k).pos0 = [5; 2 + 0.6*k];
            obs(k).vel  = [0;0];
            obs(k).R    = 0.8;
        end
    end
end

function RHS = quadrotor_rhs_collocation(states, controls, p) %#ok<INUSD>
    RHS = [ states(3,:);
            states(4,:);
            controls(1,:);
            controls(2,:) ];
end

function integrand = quadrotor_cost_integrand(states, controls, u_dot, p, w) %#ok<INUSD>
    ax  = controls(1,:); ay  = controls(2,:);
    dax = u_dot(1,:);    day = u_dot(2,:);

    switch w.objective_mode
        case 'time_energy_L2'
            integrand = w.C + w.wu(1)*(ax.^2) + w.wu(2)*(ay.^2);
        case 'min_energy_L2'
            integrand = w.wu(1)*(ax.^2) + w.wu(2)*(ay.^2);
        case 'min_jerk_L2'
            integrand = w.wudot(1)*(dax.^2) + w.wudot(2)*(day.^2);
            if any(w.wu_reg_in_jerk > 0)
                integrand = integrand + w.wu_reg_in_jerk(1)*(ax.^2) + w.wu_reg_in_jerk(2)*(ay.^2);
            end
        otherwise
            error('Unsupported objective_mode: %s', w.objective_mode);
    end
end

function c = quadrotor_ineq_constraints(states, controls, time_nodes, p, b) %#ok<INUSD>
% Hermite-UPSAMPLED inside intervals to reduce between-node tunneling
% Convention: fmincon uses c <= 0 (positive => violation)

    x  = states(1,:); y  = states(2,:);
    vx = states(3,:); vy = states(4,:);

    N = numel(time_nodes) - 1;
    nSub = max(1, p.constraint_upsample);

    tS  = []; xS  = []; yS  = [];
    vxS = []; vyS = [];

    for k = 1:N
        tk  = time_nodes(k);
        tk1 = time_nodes(k+1);
        h = (tk1 - tk);

        s = linspace(0,1,nSub+1);
        if k < N, s = s(1:end-1); end

        t_local = tk + (tk1-tk)*s;

        [x_local, vx_local] = hermite_1d(x(k),  vx(k),  x(k+1),  vx(k+1),  h, s);
        [y_local, vy_local] = hermite_1d(y(k),  vy(k),  y(k+1),  vy(k+1),  h, s);

        tS  = [tS,  t_local];  %#ok<AGROW>
        xS  = [xS,  x_local];  %#ok<AGROW>
        yS  = [yS,  y_local];  %#ok<AGROW>
        vxS = [vxS, vx_local]; %#ok<AGROW>
        vyS = [vyS, vy_local]; %#ok<AGROW>
    end

    % Speed: ||v|| <= v_max  =>  ||v|| - v_max <= 0
    speed = sqrt(vxS.^2 + vyS.^2 + 1e-12);
    c_speed = (speed - p.v_max)';

    % Map bounds (optional): keep within [xlim,ylim]
    c_map = [];
    if isfield(b,'xlim') && isfield(b,'ylim') && numel(b.xlim)==2 && numel(b.ylim)==2
        xmin = b.xlim(1); xmax = b.xlim(2);
        ymin = b.ylim(1); ymax = b.ylim(2);

        c_map = [ (xmin - xS).'; (xS - xmax).'; (ymin - yS).'; (yS - ymax).' ];
    end

    % Obstacles: dist >= (R + safety_margin)  <=>  (R+sm)^2 - dist^2 <= 0
    nObs = numel(p.obs);
    nS   = numel(tS);
    c_obs = zeros(nObs*nS, 1);

    idx = 1;
    for i = 1:nObs
        ox = p.obs(i).pos0(1) + p.obs(i).vel(1)*tS;
        oy = p.obs(i).pos0(2) + p.obs(i).vel(2)*tS;

        dist_sq = (xS-ox).^2 + (yS-oy).^2;
        min_dist = (p.obs(i).R + p.safety_margin);
        ci = (min_dist^2 - dist_sq)'; % <= 0 feasible

        c_obs(idx:idx+nS-1) = ci;
        idx = idx + nS;
    end

    % Optional control rate bound at nodes
    c_du = [];
    if isfield(b,'du_max') && ~isempty(b.du_max)
        du = diff(controls,1,2);       % 2 x N
        du_max = b.du_max(:);          % 2 x 1
        c1 = du - du_max;              % <=0
        c2 = -du - du_max;             % <=0
        c_du = [c1(:); c2(:)];
    end

    c = [c_speed; c_map(:); c_obs; c_du(:)];
    c = c(:);
end

function quadrotor_summary(solution, summary, problem)
    T_final  = solution.T_final;
    cost     = solution.cost;
    X        = solution.states;

    speed = sqrt(X(3,:).^2 + X(4,:).^2);
    max_speed = max(speed);

    fprintf('====================================================================\n');
    fprintf('QUADROTOR NAV SUMMARY | %s\n', problem.name);
    fprintf('--------------------------------------------------------------------\n');
    fprintf('T_final                : %.4f\n', T_final);
    fprintf('J_final                : %.6f\n', cost);
    fprintf('Runtime (total)        : %.2f sec\n', summary.total_time);
    fprintf('Max speed ||v||        : %.4f  (limit %.2f)\n', max_speed, problem.params.v_max);

    % Backward compatible printing (old vs new summary fields)
    if isfield(summary,'qpso_best_maxDenseIneqViol')
        fprintf('QPSO best maxDense     : %.2e\n', summary.qpso_best_maxDenseIneqViol);
        fprintf('QPSO best minDenseClear: %.3f\n', summary.qpso_best_minDenseClear);
        fprintf('QPSO best termErr      : %.2e\n', summary.qpso_best_termErr);
    else
        if isfield(summary,'qpso_best_max_viol'),  fprintf('QPSO best maxViol      : %.2e\n', summary.qpso_best_max_viol); end
        if isfield(summary,'qpso_best_min_clear'), fprintf('QPSO best minClear     : %.3f\n', summary.qpso_best_min_clear); end
        if isfield(summary,'qpso_best_term_err'),  fprintf('QPSO best termErr      : %.2e\n', summary.qpso_best_term_err); end
    end

    fprintf('PMP status             : %s\n', solution.pmp_status);
    fprintf('====================================================================\n\n');
end

function [dH_du_all, d2H_du2_all, dH_du_dt_all] = quadrotor_pmp_derivatives( ...
    states, controls, costates, p, w, time_vec) %#ok<INUSD>

    if ~any(strcmp(w.objective_mode, {'time_energy_L2','min_energy_L2'}))
        error('PMP derivatives only valid for time_energy_L2 or min_energy_L2.');
    end
    if ~isfield(w,'wu') || any(w.wu <= 0)
        error('PMP derivatives require positive L2 accel weights w.wu.');
    end

    ax = controls(1,:);
    ay = controls(2,:);
    lambda_vx = costates(3,:);
    lambda_vy = costates(4,:);

    dH_dax = 2*w.wu(1)*ax + lambda_vx;
    dH_day = 2*w.wu(2)*ay + lambda_vy;

    d2H_dax2 = 2*w.wu(1)*ones(size(ax));
    d2H_day2 = 2*w.wu(2)*ones(size(ay));

    dt = time_vec(2) - time_vec(1);
    dH_dax_dt = gradient(dH_dax, dt);
    dH_day_dt = gradient(dH_day, dt);

    dH_du_all    = [dH_dax; dH_day];
    d2H_du2_all  = [d2H_dax2; d2H_day2];
    dH_du_dt_all = [dH_dax_dt; dH_day_dt];
end

%==========================================================================
% Hermite 1D interpolation: position & velocity along an interval
%==========================================================================
function [p_s, v_s] = hermite_1d(p0, v0, p1, v1, h, s)
    s = s(:).'; % 1xM
    h00 =  2*s.^3 - 3*s.^2 + 1;
    h10 =      s.^3 - 2*s.^2 + s;
    h01 = -2*s.^3 + 3*s.^2;
    h11 =      s.^3 - s.^2;

    p_s = p0*h00 + (h*v0)*h10 + p1*h01 + (h*v1)*h11;

    dh00 = 6*s.^2 - 6*s;
    dh10 = 3*s.^2 - 4*s + 1;
    dh01 = -dh00;
    dh11 = 3*s.^2 - 2*s;

    dp_ds = p0*dh00 + (h*v0)*dh10 + p1*dh01 + (h*v1)*dh11;
    v_s = (1/h) * dp_ds;
end

function v = getfield_default(s, f, d)
    if isfield(s,f) && ~isempty(s.(f)), v = s.(f);
    else, v = d; end
end
