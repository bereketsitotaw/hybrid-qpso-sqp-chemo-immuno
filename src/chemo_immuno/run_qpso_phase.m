function [gBest, gBestCost, qpso_time, qpso] = run_qpso_phase(p, N, T_guess, qpso)
%RUN_QPSO_PHASE Stage-1 global search with Quantum PSO (Section 4).

    num_control_vars = 2 * (N + 1);
    lb_controls = [zeros(1, N + 1), zeros(1, N + 1)];
    ub_controls = [p.umax * ones(1, N + 1), p.vmax * ones(1, N + 1)];

    particles = lb_controls + (ub_controls - lb_controls) .* rand(qpso.swarmSize, num_control_vars);
    pBest = particles;
    pBestCost = inf(qpso.swarmSize, 1);
    gBest = particles(1, :);
    gBestCost = inf;

    tic;
    for iter = 1:qpso.maxIter
        if qpso.useParfor
            parfor i = 1:qpso.swarmSize
                u = particles(i, 1:N + 1);
                v = particles(i, N + 2:end);
                cost = qpso_cost_func(u, v, T_guess, p, N, qpso.rho, qpso);
                if cost < pBestCost(i)
                    pBestCost(i) = cost;
                    pBest(i, :) = particles(i, :);
                end
            end
        else
            for i = 1:qpso.swarmSize
                u = particles(i, 1:N + 1);
                v = particles(i, N + 2:end);
                cost = qpso_cost_func(u, v, T_guess, p, N, qpso.rho, qpso);
                if cost < pBestCost(i)
                    pBestCost(i) = cost;
                    pBest(i, :) = particles(i, :);
                end
            end
        end

        [min_pBestCost, idx] = min(pBestCost);
        if min_pBestCost < gBestCost
            gBestCost = min_pBestCost;
            gBest = pBest(idx, :);
        end

        if qpso.verbose
            fprintf('      QPSO Iteration: %d / %d | Best Cost Found: %.4f\n', ...
                iter, qpso.maxIter, gBestCost);
        end

        mBest = mean(pBest, 1);
        for i = 1:qpso.swarmSize
            phi = rand(1, num_control_vars);
            p_i = phi .* pBest(i, :) + (1 - phi) .* gBest;
            u_val = rand(1, num_control_vars);
            L = 2 * abs(particles(i, :) - mBest);

            noise_amplitude = qpso.noise_amplitude;
            eta = randn(1, num_control_vars);
            uniform_noise = rand(1, num_control_vars) - 0.5;
            exploration_term = noise_amplitude * eta .* uniform_noise;

            particles(i, :) = p_i + (-1)^ceil(0.5 + rand) .* L .* log(1 ./ u_val) + exploration_term;
        end
        particles = max(min(particles, ub_controls), lb_controls);
    end
    qpso_time = toc;
end
