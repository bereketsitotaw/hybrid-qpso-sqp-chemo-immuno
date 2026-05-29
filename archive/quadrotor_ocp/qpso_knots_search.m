%==========================================================================
% qpso_knots_search.m
%==========================================================================
% QPSO over low-dimensional control knots with:
%   - feasibility-first scoring using dense inequality checks
%   - continuation schedule for penalty emphasis
%   - restart logic
%   - optional adaptive knot count increase when infeasible
%==========================================================================
function [knotsBest, U_best, stats] = qpso_knots_search(problem, T_guess, tr, qpso)

    N  = problem.time.N;
    nc = problem.num_controls;

    umin = problem.bounds.control_min(:);
    umax = problem.bounds.control_max(:);

    % continuation schedule (soft -> strict)
    stageSched = qpso_penalty_schedule(qpso.nStages);

    % stats
    stats = struct();
    stats.qpso_evals = 0;
    stats.bestCost = inf;
    stats.bestMetrics = struct('maxDenseIneqViolation', inf, 'minDenseClearance', -inf, 'termErr', inf);

    % outer loop: allow knot adaptation (rarely triggers)
    nKnots = qpso.nKnots;

    while true
        num_u_vars = nc*nKnots;

        % Knot-space bounds (flat)
        lbK  = repmat(umin, nKnots, 1);
        ubK  = repmat(umax, nKnots, 1);

        fprintf('   [QPSO] nKnots=%d | swarm=%d | iters=%d | stages=%d\n', ...
            nKnots, qpso.swarmSize, qpso.maxIter, qpso.nStages);

        % init swarm uniformly
        particles = lbK.' + (ubK.'-lbK.').*rand(qpso.swarmSize, num_u_vars);

        pBest     = particles;
        pBestCost = inf(qpso.swarmSize,1);
        gBest     = particles(1,:);
        gBestCost = inf;
        gBestMetrics = stats.bestMetrics;

        for stage = 1:qpso.nStages
            wpen = stageSched(stage);
            fprintf('   QPSO stage %d/%d | A(viol)=%.1e | B(term)=%.1e\n', ...
                stage, qpso.nStages, wpen.A_viol, wpen.B_term);

            stallCount = 0;
            gBestCost_stage_start = gBestCost;

            % ramp: if no feasible found, increase A_viol gradually
            A_base = wpen.A_viol;

            for iter = 1:qpso.maxIter

                A_viol = A_base * (1 + 5*(iter/qpso.maxIter));  % ramps within stage

                cost_i = zeros(qpso.swarmSize,1);

                if qpso.useParfor
                    parfor i = 1:qpso.swarmSize
                        cost_i(i) = qpso_cost_knots(particles(i,:), T_guess, problem, tr, nKnots, A_viol, wpen.B_term);
                    end
                else
                    for i = 1:qpso.swarmSize
                        cost_i(i) = qpso_cost_knots(particles(i,:), T_guess, problem, tr, nKnots, A_viol, wpen.B_term);
                    end
                end
                stats.qpso_evals = stats.qpso_evals + qpso.swarmSize;

                improved = cost_i < pBestCost;
                pBestCost(improved) = cost_i(improved);
                pBest(improved,:)   = particles(improved,:);

                [minCost, idx] = min(pBestCost);
                if minCost < gBestCost
                    gBestCost = minCost;
                    gBest     = pBest(idx,:);
                    stallCount = 0;

                    % refresh best metrics
                    [~, m] = qpso_cost_knots(gBest, T_guess, problem, tr, nKnots, A_viol, wpen.B_term);
                    gBestMetrics = m;
                else
                    stallCount = stallCount + 1;
                end

                if mod(iter,20)==0 || iter==1 || iter==qpso.maxIter
                    fprintf('      iter %3d/%d | bestCost %.3e | maxDense %.2e | minClear %.3f | termErr %.2e\n', ...
                        iter, qpso.maxIter, gBestCost, gBestMetrics.maxDenseIneqViolation, ...
                        gBestMetrics.minDenseClearance, gBestMetrics.termErr);
                end

                % restart
                if stallCount >= qpso.stallIters
                    nR = max(1, round(qpso.restartFrac*qpso.swarmSize));
                    ridx = randperm(qpso.swarmSize, nR);
                    spread = 0.25;
                    particles(ridx,:) = gBest + spread*randn(nR, num_u_vars).*repmat((ubK-lbK).', nR, 1);
                    stallCount = 0;
                end

                % QPSO update
                beta_qpso = qpso.beta0 + (qpso.betaF - qpso.beta0)*(iter/qpso.maxIter);
                noise_amp = qpso.noise0 + (qpso.noiseF - qpso.noise0)*(iter/qpso.maxIter);
                mBest = mean(pBest,1);

                for i = 1:qpso.swarmSize
                    phi = rand(1,num_u_vars);
                    p_i = phi.*pBest(i,:) + (1-phi).*gBest;

                    u = rand(1,num_u_vars);
                    L = abs(mBest - particles(i,:));

                    signVec = 2*(rand(1,num_u_vars) > 0.5) - 1; % per-dim sign (FIX)
                    step = signVec .* beta_qpso .* L .* log(1./max(u,1e-12));
                    exploration = noise_amp * randn(1,num_u_vars);

                    particles(i,:) = p_i + step + exploration;
                end

                particles = max(min(particles, ubK.'), lbK.');
            end

            if gBestCost >= gBestCost_stage_start
                fprintf('   (note) stage %d did not improve global best.\n', stage);
            end
        end

        % If requested, adapt knots when still clearly infeasible
        if qpso.adaptKnots && (gBestMetrics.minDenseClearance < -qpso.adaptClearEps)
            fprintf('   [QPSO] Still infeasible (minClear=%.3e). Increasing knots %d -> %d and retrying.\n', ...
                gBestMetrics.minDenseClearance, nKnots, nKnots + qpso.adaptKnotsStep);
            nKnots = nKnots + qpso.adaptKnotsStep;

            % brief exit if knots become huge
            if nKnots > 60
                warning('[QPSO] nKnots exceeded 60; stopping adaptation.');
                break;
            end
        else
            % done
            knotsBest = gBest(:);
            U_best = knots_to_controls(knotsBest, T_guess, problem, nKnots); %#ok<NASGU>
            stats.bestCost = gBestCost;
            stats.bestMetrics = gBestMetrics;
            break;
        end
    end

    % final outputs
    knotsBest = gBest(:);
    U_best = knots_to_controls(knotsBest, T_guess, problem, nKnots);
    stats.bestCost = gBestCost;
    stats.bestMetrics = gBestMetrics;
    stats.nKnotsFinal = nKnots;
end

%-------------------- schedule --------------------
function sched = qpso_penalty_schedule(nStages)
    if nStages < 1, nStages = 1; end
    sched = repmat(struct('A_viol',1e4,'B_term',1e3), 1, nStages);

    if nStages == 1
        sched(1) = struct('A_viol',1e6,'B_term',1e4);
        return;
    end

    sched(1) = struct('A_viol',1e3,'B_term',1e2);
    if nStages >= 2
        sched(2) = struct('A_viol',1e5,'B_term',1e3);
    end
    for s = 3:nStages
        sched(s) = struct('A_viol',1e7,'B_term',1e4);
    end
end

%-------------------- cost + metrics --------------------
function [cost, m] = qpso_cost_knots(knots_flat, T, problem, tr, nKnots, A_viol, B_term)

    nc = problem.num_controls;

    U = knots_to_controls(knots_flat(:), T, problem, nKnots);
    x0 = problem.bounds.initial_states(:);

    X = tr.rollout(x0, U, T);

    % guardrails: blow-up => huge cost
    if any(~isfinite(X(:))) || any(isnan(X(:)))
        cost = 1e12;
        m = struct('maxDenseIneqViolation',1e6,'minDenseClearance',-1e6,'termErr',1e6);
        return;
    end

    mt = tr.metrics(X, U, T);

    viol = max(0, mt.maxDenseIneqViolation);
    cost = A_viol*viol + B_term*mt.termErr + mt.J;

    if ~isfinite(cost) || isnan(cost)
        cost = 1e12;
    end

    m = struct('maxDenseIneqViolation', mt.maxDenseIneqViolation, ...
               'minDenseClearance', mt.minDenseClearance, ...
               'termErr', mt.termErr);
end

%-------------------- knots -> node controls --------------------
function U = knots_to_controls(knots_flat, T, problem, nKnots)
    N  = problem.time.N;
    nc = problem.num_controls;

    knots = reshape(knots_flat, nc, nKnots);
    tk = linspace(0, T, nKnots);
    tn = linspace(0, T, N+1);

    U = zeros(nc, N+1);
    for i = 1:nc
        U(i,:) = interp1(tk, knots(i,:), tn, 'pchip', 'extrap');
    end

    umin = problem.bounds.control_min(:);
    umax = problem.bounds.control_max(:);
    for i = 1:nc
        U(i,:) = min(max(U(i,:), umin(i)), umax(i));
    end
end
