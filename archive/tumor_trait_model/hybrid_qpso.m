%==========================================================================
% Optimal Control of Tumor Population Dynamics Using Hybrid QPSO
%==========================================================================
%% 1. Configuration & Setup
%--------------------------------------------------------------------------
clear;
clc;
close all;
tic; % Start timer
rng(42); % Set seed for reproducibility

fprintf('--- QPSO Optimization (Nominal Benchmark) ---\n');

% Load shared configuration
cfg = config();

% --- Simulation & Model Parameters ---
p.T = cfg.model.T;
p.n = cfg.model.n_steps;
p.m = cfg.model.m;
p.time = linspace(0, p.T, p.n + 1);
p.x_trait = cfg.model.x;
p.e = cfg.model.e;

% Matrix forms for dynamics
p.R = cfg.model.R;
p.Phi1 = cfg.model.Phi; % Using shared susceptibility model
p.Phi2 = zeros(p.m, p.m); % Drug 2 disabled for comparison
p.M = cfg.model.M;

% --- Mutation & Interaction Parameters ---
% Force NO MUTATION for fair open-loop comparison with PMP
p.theta = 0.0;
p.sigma = 5;
p.interaction = 0.0;
p.Phi12 = zeros(p.m, p.m); % No interaction

% --- Cost Function Weights ---
% Exact match to PMP's vectorized scaling: alpha_vec = alpha * m * e
p.alpha_vec = cfg.cost.alpha * p.m * p.e; 
p.beta_vec  = cfg.cost.beta * p.m * p.e;
p.gamma1    = cfg.cost.gamma;
p.gamma2    = 0;

% --- Control Constraints ---
p.Umax1 = cfg.cost.Umax;
p.Umax2 = 0;
num_drugs = 1;

% --- QPSO & fmincon Settings ---
opt.qpso.swarmSize = 50; 
opt.qpso.maxIter = 100;
opt.fmincon = optimoptions('fmincon', 'Algorithm', 'sqp', ...
    'Display', 'iter', ...
    'MaxIterations', 50000, ...
    'MaxFunctionEvaluations', 1e6);
    
%% 2. Derived Parameters & Initial Conditions
%--------------------------------------------------------------------------
n0 = cfg.model.N0;
lb1 = zeros(1, p.n + 1);
ub1 = p.Umax1 * ones(1, p.n + 1);
lb2 = zeros(1, p.n + 1);
ub2 = p.Umax2 * ones(1, p.n + 1); % This will be all zeros for the single-drug case
lb_full = [lb1, lb2];
ub_full = [ub1, ub2];

%% 3. Optimization Phase
%--------------------------------------------------------------------------
% (No changes in this section)
fprintf('Starting QPSO optimization phase...\n');
particles1 = lb1 + (ub1 - lb1) .* rand(opt.qpso.swarmSize, p.n + 1);
particles2 = lb2 + (ub2 - lb2) .* rand(opt.qpso.swarmSize, p.n + 1);
pBest1 = particles1;
pBest2 = particles2;
pBestCost = inf(opt.qpso.swarmSize, 1);
gBest1 = particles1(1, :);
gBest2 = particles2(1, :);
gBestCost = inf;
for iter = 1:opt.qpso.maxIter
    parfor i = 1:opt.qpso.swarmSize
        [~, ~, totalCost] = computeCost(particles1(i, :), particles2(i, :), p, n0);
        if totalCost < pBestCost(i)
            pBestCost(i) = totalCost;
            pBest1(i, :) = particles1(i, :);
            pBest2(i, :) = particles2(i, :);
        end
    end
    [min_pBestCost, idx] = min(pBestCost);
    if min_pBestCost < gBestCost
        gBestCost = min_pBestCost;
        gBest1 = pBest1(idx, :);
        gBest2 = pBest2(idx, :);
    end
    fprintf('QPSO Iteration: %d / %d | Best Cost: %.2f\n', iter, opt.qpso.maxIter, gBestCost);
    p_mean1 = mean(pBest1, 1);
    p_mean2 = mean(pBest2, 1);
    ksi1 = rand(opt.qpso.swarmSize, p.n + 1);
    ksi2 = rand(opt.qpso.swarmSize, p.n + 1);
    eta1 = randn(opt.qpso.swarmSize, p.n + 1);
    eta2 = randn(opt.qpso.swarmSize, p.n + 1);
    particles1 = p_mean1 + ksi1 .* abs(gBest1 - particles1) .* log(1 ./ rand(size(particles1))) + eta1 .* (rand(size(particles1)) - 0.5);
    particles2 = p_mean2 + ksi2 .* abs(gBest2 - particles2) .* log(1 ./ rand(size(particles2))) + eta2 .* (rand(size(particles2)) - 0.5);
    particles1 = min(max(particles1, lb1), ub1);
    particles2 = min(max(particles2, lb2), ub2);
end
fprintf('QPSO phase complete. Best cost found: %.2f\n\n', gBestCost);
fprintf('Starting fmincon refinement phase...\n');
costFunc = @(u) computeCostWrapper(u, p, n0);
initialGuess = [gBest1, gBest2];
[optimal_controls_vec, finalCost] = fmincon(costFunc, initialGuess, [], [], [], [], lb_full, ub_full, [], opt.fmincon);
fprintf('fmincon refinement complete. Final cost: %.2f\n', finalCost);

%% 4. Post-Processing & Plotting
%--------------------------------------------------------------------------
% (Minor change in this section for plotting)
fprintf('\nComputing final trajectories and generating plots...\n');
u1_opt_vec = optimal_controls_vec(1 : p.n+1);
u2_opt_vec = optimal_controls_vec(p.n+2 : end);
u1_func = @(t) interp1(p.time, u1_opt_vec, t, 'pchip');
u2_func = @(t) interp1(p.time, u2_opt_vec, t, 'pchip');
N_final_traj = rk4Method(u1_func, u2_func, p, n0);
fprintf('\nSaving QPSO results to qpso_results.mat...\n');
qpso_results = struct();
qpso_results.time = p.time;
qpso_results.u1 = u1_opt_vec;
qpso_results.u2 = u2_opt_vec;
qpso_results.N = N_final_traj;
save('qpso_results.mat', 'qpso_results');
fprintf('Results saved successfully.\n');
figure('Name', 'Optimal Control Results', 'Position', [100, 100, 1200, 800]);
sgtitle('Optimal Control Results (Nominal QPSO Benchmark)', 'FontSize', 16, 'FontWeight', 'bold');
subplot(2, 2, 1);
hold on;
plot(p.time, u1_opt_vec, 'b-', 'LineWidth', 2, 'DisplayName', 'Drug 1 (u_1)');
% --- MODIFIED: Only plot drug 2 if it's being used ---
if num_drugs > 1
    plot(p.time, u2_opt_vec, 'r--', 'LineWidth', 2, 'DisplayName', 'Drug 2 (u_2)');
end
hold off;
title('Optimal Control Inputs'); xlabel('Time (years)'); ylabel('Dosage');
grid on; legend show; ylim([-0.1, p.Umax1 + 0.2]);

% (The rest of the plotting is unchanged)
subplot(2, 2, 2);
plot(p.time, mean(N_final_traj, 1), 'k-', 'LineWidth', 2);
title('Average Tumor Population'); xlabel('Time (years)'); ylabel('Avg. Population');
grid on;
subplot(2, 2, 3);
plot(p.x_trait, N_final_traj(:, 1), 'k--', 'LineWidth', 2, 'DisplayName', 'Initial');
hold on;
plot(p.x_trait, N_final_traj(:, end), 'm-', 'LineWidth', 2, 'DisplayName', 'Final');
hold off;
title('Tumor Distribution Across Traits'); xlabel('Resistance Trait (x)'); ylabel('Population');
grid on; legend show;
subplot(2, 2, 4);
surf(p.time, p.x_trait, N_final_traj, 'EdgeColor', 'none');
title('Tumor Population Evolution'); xlabel('Time (years)'); ylabel('Resistance Trait (x)');
zlabel('Population'); view(30, 45); colorbar;

%% 5. Display Performance Summary
%--------------------------------------------------------------------------
% (No changes in this section)
% Use shared compute_metrics script ensuring numerical equality with PMP evaluation
metrics = compute_metrics(p.time, N_final_traj, u1_opt_vec, cfg);

fprintf('\n\n================== PERFORMANCE SUMMARY (QPSO) ==================\n');
fprintf(' Agent Type:         QPSO\n');
fprintf('--------------------------------------------------------------------\n');
fprintf(' Final Avg Tumor:      %.3f\n', metrics.FinalTumor / p.m);
fprintf(' Tumor Reduction:      %.2f %%\n', metrics.TumorReduction);
fprintf('--------------------------------------------------------------------\n');
fprintf(' Total Drug Used:      %.3f\n', metrics.TotalDrug);
fprintf('--------------------------------------------------------------------\n');
fprintf(' Total Cost (J):       %.2f\n', metrics.TotalCost);
fprintf('   - Terminal Cost:    %.2f\n', metrics.TerminalCost);
fprintf('   - Running Tumor:    %.2f\n', metrics.RunningTumorCost);
fprintf('   - Drug Cost:        %.2f\n', metrics.DrugCost);
fprintf('====================================================================\n');
total_time = toc;
fprintf('\nTotal execution time: %.2f seconds.\n', total_time);

%% 6. Local Function Definitions
%--------------------------------------------------------------------------
% (No changes in this section)
function totalCost = computeCostWrapper(u_vec, p, n0)
    u1_vec = u_vec(1 : p.n+1);
    u2_vec = u_vec(p.n+2 : end);
    [~, ~, totalCost] = computeCost(u1_vec, u2_vec, p, n0);
end
function [instantaneousCost, terminalCost, totalCost] = computeCost(u1_vec, u2_vec, p, n0)
    u1 = @(t) interp1(p.time, u1_vec, t, 'pchip');
    u2 = @(t) interp1(p.time, u2_vec, t, 'pchip');
    N = rk4Method(u1, u2, p, n0);
    
    % Vectorized exact match to PMP
    running_tumor_cost = p.beta_vec' * N;
    integrand = running_tumor_cost + p.gamma1 * u1(p.time) + p.gamma2 * u2(p.time);
    
    instantaneousCost = trapz(p.time, integrand);
    terminalCost = p.alpha_vec' * N(:, end);
    
    % Smoothness regularization to suppress end-of-horizon control wiggles
    dt = p.time(2) - p.time(1);
    lambda_smooth = 1e2;   % initial tuning value
    smoothnessCost = lambda_smooth * sum((diff(u1_vec)/dt).^2) * dt;

    totalCost = instantaneousCost + terminalCost + smoothnessCost;
end
function N = rk4Method(u1, u2, p, n0)
    dt = p.time(2) - p.time(1);
    N = zeros(p.m, p.n + 1);
    N(:, 1) = n0;
    
    dyn_params.e = p.e;
    dyn_params.R = p.R;
    dyn_params.Phi = p.Phi1;
    dyn_params.M = p.M;

    for k = 1:p.n
        t = p.time(k);
        N_current = N(:, k);
        u_val = max(0, min(u1(t), p.Umax1));
        u_val_mid = max(0, min(u1(t + dt/2), p.Umax1));
        u_val_next = max(0, min(u1(t + dt), p.Umax1));
        
        k1 = tumor_dynamics(N_current, u_val, dyn_params);
        k2 = tumor_dynamics(N_current + dt/2*k1, u_val_mid, dyn_params);
        k3 = tumor_dynamics(N_current + dt/2*k2, u_val_mid, dyn_params);
        k4 = tumor_dynamics(N_current + dt*k3, u_val_next, dyn_params);
        
        N(:, k + 1) = N_current + (dt/6) * (k1 + 2*k2 + 2*k3 + k4);
    end
end
% analyze_performance helper deleted in favor of compute_metrics.m