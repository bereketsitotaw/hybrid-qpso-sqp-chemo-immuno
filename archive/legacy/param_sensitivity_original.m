%==========================================================================
% MATLAB Code for Sensitivity Analysis & Pareto Front Generation
% Method: Hybrid QPSO + fmincon (Direct Method)
% Final, Merged & Enhanced Version for Publication
%==========================================================================
%% 1. Initialization and Setup
clear; clc; close all;
rng(1); % for reproducibility

% --- Discretization and Time ---
N = 50; % Reduced grid points for faster analysis
T_guess = 4.0; % Initial guess for final time T

% --- Model Parameters (Structure 'p') ---
p.xi = 0.5618;    p.x_inf = 780;    p.theta = 1.0;    p.alpha = 1.0;
p.a = 0.00726;    p.b = 0.00264;    p.gamma = 0.1181; p.delta = 0.3745;
p.kappa = 2.0;    p.nu = 1.0;       p.sigma = 0.1;    p.psi = 2.0;
p.mu = 0.45;

% --- Boundary Conditions and Control Limits ---
p.x0 = 650;       p.y0 = 0.1;       p.z0 = 400;       p.xf = 1.0;
p.umax = 1.0;     p.vmax = 1.0;

% --- QPSO Parameters ---
qpso.swarmSize = 50; % Reduced for faster analysis runs
qpso.maxIter = 50;   % Reduced for faster analysis runs

% --- Define weight ranges for sensitivity analysis ---
l1_range = [0, 0.1, 1, 10];
l2_range = [0, 0.1, 1, 10];
l3_range = [0, 0.1, 1, 10];
A_range  = [0, 0.1, 0.3, 1, 3, 10, 30, 50]; % Expanded for Pareto
B_range  = [0, 0.3, 1, 10, 50];
C_range  = [0, 1, 5, 10, 50];

% --- Initialize storage for results ---
results_l1 = cell(1, length(l1_range));
results_l2 = cell(1, length(l2_range));
results_l3 = cell(1, length(l3_range));
results_A  = cell(1, length(A_range));
results_B  = cell(1, length(B_range));
results_C  = cell(1, length(C_range));

%% 2. Sensitivity Analysis Loops
% --- Define the Baseline Weights ---
baseline_p = p;
baseline_p.l1 = 0;   baseline_p.l2 = 0;   baseline_p.l3 = 0; 
baseline_p.A = 0;    baseline_p.B = 0;    baseline_p.C = 1.0;

% --- Run sensitivity sweeps by varying one parameter at a time ---
fprintf('Running sensitivity analysis...\n');
all_ranges = {l1_range, l2_range, l3_range, A_range, B_range, C_range};
all_results = {results_l1, results_l2, results_l3, results_A, results_B, results_C};
param_fields = {'l1', 'l2', 'l3', 'A', 'B', 'C'};

for k = 1:length(param_fields)
    fprintf('Sweeping parameter: %s\n', param_fields{k});
    current_p = baseline_p;
    current_range = all_ranges{k};
    temp_results = cell(1, length(current_range));
    
    parfor i = 1:length(current_range)
        p_loop = current_p;
        p_loop.(param_fields{k}) = current_range(i);
        [z_opt, J_opt, time_opt] = run_optimization(p_loop, N, T_guess, qpso);
        
        temp_results{i} = struct(...
            'J_opt', J_opt, ...
            'T_opt', z_opt(end), ...
            'u_opt', z_opt(3*(N+1)+1 : 4*(N+1)), ...
            'v_opt', z_opt(4*(N+1)+1 : 5*(N+1)), ...
            'time_opt', time_opt, ...
            'param_val', current_range(i) ...
        );
        
        fprintf('  %s = %.2f, J_opt = %.4f, T_opt = %.2f\n', param_fields{k}, current_range(i), J_opt, z_opt(end));
    end
    all_results{k} = temp_results;
end

results_l1 = all_results{1}; results_l2 = all_results{2}; results_l3 = all_results{3};
results_A = all_results{4}; results_B = all_results{5}; results_C = all_results{6};

%% 3. Data Preparation for Advanced Visualization
fprintf('Preparing data for advanced plots...\n');

% Data for Chemo Panel (vs. A)
total_u_dose_A = cellfun(@(r) trapz(r.time_opt, r.u_opt), results_A);
total_v_dose_A = cellfun(@(r) trapz(r.time_opt, r.v_opt), results_A);
T_opts_A = cellfun(@(r) r.T_opt, results_A);
J_opts_A = cellfun(@(r) r.J_opt, results_A);

% Data for Immuno Panel (vs. B)
total_u_dose_B = cellfun(@(r) trapz(r.time_opt, r.u_opt), results_B);
total_v_dose_B = cellfun(@(r) trapz(r.time_opt, r.v_opt), results_B);
T_opts_B = cellfun(@(r) r.T_opt, results_B);
J_opts_B = cellfun(@(r) r.J_opt, results_B);

% Data for Tornado Plot
param_names = {'$l_1$ (Tumor)', '$l_2$ (Immune)', '$l_3$ (Antigen)', '$A$ (Chemo)', '$B$ (Immuno)', '$C$ (Time)'};
sensitivities = zeros(length(param_names), 1);
results_all_struct = {results_l1, results_l2, results_l3, results_A, results_B, results_C};
for i = 1:length(param_names)
    J_vals = cellfun(@(r) r.J_opt, results_all_struct{i});
    J_baseline = J_vals(1);
    J_max = max(J_vals);
    if abs(J_baseline) < 1e-9; sensitivities(i) = 100 * J_max; else; sensitivities(i) = (J_max - J_baseline) / J_baseline * 100; end
end
[sorted_sens, sort_idx] = sort(sensitivities, 'descend');
sorted_names = param_names(sort_idx);

%% 4. Advanced Visualization & Publication-Quality Export
fprintf('Generating publication-quality plots...\n');

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');
set(groot, 'defaultAxesFontSize', 14);
plotLineWidth = 2.5; markerSize = 8;

% --- PLOT 1: Chemo Sensitivity Panel (vs. A) ---
fig1 = figure('Name', 'Chemo Sensitivity Panel', 'Position', [100, 100, 1000, 800]);
sgtitle('Sensitivity Analysis with Respect to Chemotherapy Cost ($A$)', 'FontSize', 16, 'FontWeight', 'bold');
subplot(2, 2, 1); plot(A_range, total_u_dose_A, 'o-', 'LineWidth', plotLineWidth, 'MarkerSize', markerSize, 'Color', '#0072BD'); title('Chemotherapy Usage'); ylabel('Total Chemo Dose $\int u(t)dt$'); xlabel('Chemotherapy Cost Weight ($A$)'); grid on; ylim([0, inf]);
subplot(2, 2, 2); plot(A_range, T_opts_A, '^-k', 'LineWidth', plotLineWidth, 'MarkerSize', markerSize); title('Optimal Treatment Time'); ylabel('Optimal Final Time $T_{opt}$'); xlabel('Chemotherapy Cost Weight ($A$)'); grid on; ylim([0, inf]);
subplot(2, 2, 3); plot(A_range, total_v_dose_A, 's--', 'LineWidth', plotLineWidth, 'MarkerSize', markerSize, 'Color', '#D95319'); title('Immunotherapy Usage'); ylabel('Total Immuno Dose $\int v(t)dt$'); xlabel('Chemotherapy Cost Weight ($A$)'); grid on; ylim([0, inf]);
subplot(2, 2, 4); plot(A_range, J_opts_A, 'd-r', 'LineWidth', plotLineWidth, 'MarkerSize', markerSize); title('Overall Optimal Cost'); ylabel('Optimal Cost ($J_{opt}$)'); xlabel('Chemotherapy Cost Weight ($A$)'); grid on; set(gca, 'YScale', 'log');

% --- NEW PLOT 2: Immuno Sensitivity Panel (vs. B) ---
fig2 = figure('Name', 'Immuno Sensitivity Panel', 'Position', [120, 120, 1000, 800]);
sgtitle('Sensitivity Analysis with Respect to Immunotherapy Cost ($B$)', 'FontSize', 16, 'FontWeight', 'bold');
subplot(2, 2, 1); plot(B_range, total_v_dose_B, 'o-', 'LineWidth', plotLineWidth, 'MarkerSize', markerSize, 'Color', '#77AC30'); title('Immunotherapy Usage'); ylabel('Total Immuno Dose $\int v(t)dt$'); xlabel('Immunotherapy Cost Weight ($B$)'); grid on; ylim([0, inf]);
subplot(2, 2, 2); plot(B_range, T_opts_B, '^-k', 'LineWidth', plotLineWidth, 'MarkerSize', markerSize); title('Optimal Treatment Time'); ylabel('Optimal Final Time $T_{opt}$'); xlabel('Immunotherapy Cost Weight ($B$)'); grid on; ylim([0, inf]);
subplot(2, 2, 3); plot(B_range, total_u_dose_B, 's--', 'LineWidth', plotLineWidth, 'MarkerSize', markerSize, 'Color', '#EDB120'); title('Chemotherapy Usage'); ylabel('Total Chemo Dose $\int u(t)dt$'); xlabel('Immunotherapy Cost Weight ($B$)'); grid on; ylim([0, inf]);
subplot(2, 2, 4); plot(B_range, J_opts_B, 'd-r', 'LineWidth', plotLineWidth, 'MarkerSize', markerSize); title('Overall Optimal Cost'); ylabel('Optimal Cost ($J_{opt}$)'); xlabel('Immunotherapy Cost Weight ($B$)'); grid on; set(gca, 'YScale', 'log');

% --- PLOT 3: Strategy vs. Cost Comparison ---
fig3 = figure('Name', 'Strategy vs. Cost Comparison', 'Position', [150, 150, 1200, 500]);
sgtitle('Evolution of Optimal Strategy with Increasing Drug Cost', 'FontSize', 16, 'FontWeight', 'bold');
subplot(1, 2, 1); hold on;
colors_u = parula(length(A_range));
for i = 1:length(A_range); plot(results_A{i}.time_opt, results_A{i}.u_opt, 'LineWidth', plotLineWidth, 'Color', colors_u(i,:)); end
title('Optimal Chemotherapy Dose $u(t)$'); xlabel('Time'); ylabel('Control Input $u$'); grid on;xlim([0, 4.5]); ylim([-0.1, p.umax + 0.1]); hold off;
colormap(gca, colors_u); cb_u = colorbar; cb_u.Label.String = 'Chemotherapy Cost Weight ($A$)'; cb_u.Ticks = linspace(0, 1, length(A_range)); cb_u.TickLabels = sprintfc('%g', A_range);
subplot(1, 2, 2); hold on;
colors_v = parula(length(B_range));
for i = 1:length(B_range); plot(results_B{i}.time_opt, results_B{i}.v_opt, 'LineWidth', plotLineWidth, 'Color', colors_v(i,:)); end
title('Optimal Immunotherapy Dose $v(t)$'); xlabel('Time'); ylabel('Control Input $v$'); grid on; xlim([0, 4.5]); ylim([-0.1, p.vmax + 0.1]); hold off;
colormap(gca, colors_v); cb_v = colorbar; cb_v.Label.String = 'Immunotherapy Cost Weight ($B$)'; cb_v.Ticks = linspace(0, 1, length(B_range)); cb_v.TickLabels = sprintfc('%g', B_range);

% --- PLOT 4: Tornado Plot ---
fig4 = figure('Name', 'Tornado Plot', 'Position', [200, 200, 900, 600]);
b = barh(sorted_sens); b.FaceColor = '#0072BD';
set(gca, 'YDir', 'reverse'); yticks(1:length(sorted_names)); yticklabels(sorted_names);
xlabel('Percent Change in Optimal Cost ($J_{opt}$) from Baseline'); title('Overall Sensitivity of Cost Function Weights');
grid on; set(gca, 'XScale', 'log'); ax = gca; ax.XGrid = 'on'; ax.YGrid = 'off';
xtips = b.XEndPoints; ytips = b.YData; labels = sprintfc('%.1f\\%%', b.XData);
text(xtips + 0.02*max(xtips), ytips, labels, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left', 'FontSize', 12);

% --- NEW PLOT 5: Pareto Front ---
fig5 = figure('Name', 'Pareto Front', 'Position', [250, 250, 900, 700]);
[sorted_doses, sort_idx] = sort(total_u_dose_A);
sorted_times = T_opts_A(sort_idx);
sorted_A = A_range(sort_idx);
plot(sorted_doses, sorted_times, 'o-', 'Color', '#D95319', 'LineWidth', 2.5, 'MarkerSize', 10, 'MarkerFaceColor', 'w');
grid on;
title('Pareto Front: Treatment Time vs. Chemotherapy Toxicity');
xlabel('Total Chemotherapy Dose ($\int u(t)dt$)');
ylabel('Optimal Treatment Time ($T_{opt}$ years)');
for i = 1:length(sorted_A)
    text(sorted_doses(i), sorted_times(i), sprintf('  $A = %g$', sorted_A(i)), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'FontSize', 12);
end

%% 5. Export Figures for Publication
% fprintf('Exporting figures to vector format...\n');
% if ~exist('figures', 'dir'); mkdir('figures'); end
% exportgraphics(fig1, 'figures/sensitivity_panel_A.eps', 'ContentType', 'vector');
% exportgraphics(fig2, 'figures/sensitivity_panel_B.eps', 'ContentType', 'vector');
% exportgraphics(fig3, 'figures/strategy_comparison.eps', 'ContentType', 'vector');
% exportgraphics(fig4, 'figures/tornado_plot.eps', 'ContentType', 'vector');
% exportgraphics(fig5, 'figures/pareto_front.eps', 'ContentType', 'vector');
% fprintf('Done.\n');

%% 6. Local Function: run_optimization
function [z_opt, J_opt, time_opt] = run_optimization(p, N, T_guess, qpso)
    num_control_vars = 2 * (N + 1);
    lb_controls = [zeros(1, N + 1), zeros(1, N + 1)];
    ub_controls = [p.umax * ones(1, N + 1), p.vmax * ones(1, N + 1)];
    particles = lb_controls + (ub_controls - lb_controls) .* rand(qpso.swarmSize, num_control_vars);
    pBest = particles; pBestCost = inf(qpso.swarmSize, 1);
    gBest = particles(1, :); gBestCost = inf;
    for iter = 1:qpso.maxIter
        temp_costs = zeros(qpso.swarmSize, 1);
        temp_pBest = pBest;
        temp_particles = particles;
        for i = 1:qpso.swarmSize
            u = temp_particles(i, 1:N+1); v = temp_particles(i, N+2:end);
            cost = qpso_cost_func(u, v, T_guess, p, N);
            temp_costs(i) = cost;
            if cost < pBestCost(i); temp_pBest(i, :) = temp_particles(i, :); end
        end
        pBestCost = min(pBestCost, temp_costs);
        pBest = temp_pBest;
        [min_pBestCost, idx] = min(pBestCost);
        if min_pBestCost < gBestCost; gBestCost = min_pBestCost; gBest = pBest(idx, :); end
        mBest = mean(pBest, 1);
        for i = 1:qpso.swarmSize
            phi = rand(1, num_control_vars);
            p_i = phi .* pBest(i, :) + (1 - phi) .* gBest;
            L = 2 * abs(particles(i, :) - mBest);
            noise_amplitude = 0.1;
            exploration_term = noise_amplitude * randn(1, num_control_vars) .* (rand(1, num_control_vars) - 0.5);
            particles(i, :) = p_i + (-1)^ceil(0.5+rand) .* L .* log(1 ./ rand(1, num_control_vars)) + exploration_term;
        end
        particles = max(min(particles, ub_controls), lb_controls);
    end
    u_guess_smooth = smoothdata(gBest(1:N+1), 'gaussian', 10);
    v_guess_smooth = smoothdata(gBest(N+2:end), 'gaussian', 10);
    x_guess = linspace(p.x0, p.xf, N + 1);
    y_guess = p.y0 * ones(1, N + 1);
    z_guess = p.z0 * ones(1, N + 1);
    z0 = [x_guess, y_guess, z_guess, u_guess_smooth, v_guess_smooth, T_guess]';
    lb = [zeros(1, 3*(N+1)), lb_controls, 0.1]';
    ub = [p.x_inf*ones(1,N+1), inf*ones(1,2*(N+1)), ub_controls, 20]';
    options = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'off', 'MaxIterations', 5000, 'MaxFunctionEvaluations', 1e6, 'OptimalityTolerance', 1e-6, 'ConstraintTolerance', 1e-6);
    [z_opt, J_opt] = fmincon(@(z) objective_func(z, p, N), z0, [], [], [], [], lb, ub, @(z) nonlinear_constraints(z, p, N), options);
    time_opt = linspace(0, z_opt(end), N + 1);
end
%% 7. Local Function Definitions
function J = objective_func(z, p, N)
    idx = 1:(N+1);
    x = z(idx); y = z(idx + (N+1)); z_vals = z(idx + 2*(N+1));
    u = z(idx + 3*(N+1)); v = z(idx + 4*(N+1)); T = z(end);
    time = linspace(0, T, N + 1);
    integrand = p.l1*x' + p.l2*y' + p.l3*z_vals' + p.A*u' + p.B*v' + p.C;
    J = trapz(time, integrand);
end
function [c, ceq] = nonlinear_constraints(z, p, N)
    idx = 1:(N+1);
    x = z(idx); y = z(idx + (N+1)); z_vals = z(idx + 2*(N+1));
    u = z(idx + 3*(N+1)); v = z(idx + 4*(N+1)); T = z(end);
    h = T / N;
    fx = @(x,y,u)      p.xi*x.*(1-x/p.x_inf) - p.theta*x.*y - p.alpha*x.*(u./(1+u));
    fy = @(x,y,z,u,v)  p.a*(1-p.b*x).*y.*z + p.gamma - p.delta*y - p.kappa*y.*(u./(1+u)) + p.nu*y.*(v./(1+v));
    fz = @(x,z,u)      p.sigma*x + p.psi*x.*(u./(1+u)) - p.mu*z;
    ceq_x = x(2:end) - x(1:end-1) - (h/2) * (fx(x(1:end-1), y(1:end-1), u(1:end-1)) + fx(x(2:end), y(2:end), u(2:end)));
    ceq_y = y(2:end) - y(1:end-1) - (h/2) * (fy(x(1:end-1), y(1:end-1), z_vals(1:end-1), u(1:end-1), v(1:end-1)) + fy(x(2:end), y(2:end), z_vals(2:end), u(2:end), v(2:end)));
    ceq_z = z_vals(2:end) - z_vals(1:end-1) - (h/2) * (fz(x(1:end-1), z_vals(1:end-1), u(1:end-1)) + fz(x(2:end), z_vals(2:end), u(2:end)));
    ceq_boundary = [x(1) - p.x0; y(1) - p.y0; z_vals(1) - p.z0; x(end) - p.xf];
    ceq = [ceq_x; ceq_y; ceq_z; ceq_boundary];
    c = [];
end
function cost = qpso_cost_func(u, v, T, p, N)
    time = linspace(0, T, N + 1);
    w0 = [p.x0; p.y0; p.z0];
    options = odeset('RelTol', 1e-3, 'AbsTol', 1e-3);
    try
        [~, w_traj] = ode45(@(t,w) dynamics(t,w,u,v,p,time), time, w0, options);
        x = w_traj(:,1)'; y = w_traj(:,2)'; z_vals = w_traj(:,3)';
        integrand = p.l1*x + p.l2*y + p.l3*z_vals + p.A*u + p.B*v + p.C;
        cost = trapz(time, integrand);
        if isempty(x) || abs(x(end) - p.xf) > 5.0; cost = inf; end
    catch
        cost = inf;
    end
end
function dwdt = dynamics(t, w, u_vec, v_vec, p, time_vec)
    x = w(1); y = w(2); z = w(3);
    u = interp1(time_vec, u_vec, t, 'linear', 'extrap');
    v = interp1(time_vec, v_vec, t, 'linear', 'extrap');
    dxdt = p.xi*x*(1-x/p.x_inf) - p.theta*x*y - p.alpha*x*(u/(1+u));
    dydt = p.a*(1-p.b*x)*y*z + p.gamma - p.delta*y - p.kappa*y.*(u/(1+u)) + p.nu*y.*(v/(1+v));
    dzdt = p.sigma*x + p.psi*x.*(u./(1+u)) - p.mu*z;
    dwdt = [dxdt; dydt; dzdt];
end