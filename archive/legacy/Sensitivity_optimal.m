%==========================================================================
% MATLAB Code for Sensitivity Analysis w.r.t. Relative Drug Cost (A/B)
%==========================================================================
%% 1. Initialization and Setup
clear; clc; close all;
rng(1); % for reproducibility

fprintf('--- SENSITIVITY ANALYSIS: RELATIVE DRUG COST (A/B) ---\n\n');

% --- Define the range of Cost Ratios to test ---
% We fix B=0.5 and vary A to change the ratio
cost_ratio_range = [0.1, 0.5, 1, 2, 5]; % Ratio A/B
B_fixed = 0.5;
A_range = B_fixed * cost_ratio_range;

% --- Initialize arrays to store results from each run ---
results.T_opt = zeros(size(A_range));
results.J_opt = zeros(size(A_range));
results.total_chemo = zeros(size(A_range));
results.total_immuno = zeros(size(A_range));

% --- Fixed Model Parameters ---
p_base.xi = 0.5618;    p_base.x_inf = 780;      p_base.theta = 1.0;
p_base.alpha = 1.0;    p_base.a = 0.00726;    p_base.b = 0.00264;
p_base.gamma = 0.1181; p_base.delta = 0.3745;  p_base.kappa = 2.0;
p_base.nu = 1.0;       p_base.sigma = 0.1;      p_base.psi = 2.0;
p_base.mu = 0.45;

% --- Cost Function Weights for "Balanced Protocol" (Case 2) ---
p_base.l1 = 1.0; p_base.l2 = 0; p_base.l3 = 0;
p_base.B = B_fixed;  p_base.C = 1.0;

% --- Boundary Conditions and Control Limits ---
p_base.x0 = 650;      p_base.y0 = 0.1;      p_base.z0 = 400;
p_base.xf = 1.0;      p_base.umax = 1.0;    p_base.vmax = 1.0;

% --- Discretization and QPSO Parameters ---
N = 100; T_guess = 5.0;
qpso.swarmSize = 100; qpso.maxIter = 150;

%% 2. Main Loop for Sensitivity Analysis
for i = 1:length(A_range)
    
    p = p_base; % Reset to base parameters
    p.A = A_range(i); % Set the chemo cost 'A' for this iteration
    
    fprintf('\n=====================================================\n');
    fprintf('RUN %d/%d: Solving for A/B Ratio = %.2f (A=%.2f, B=%.2f)\n', i, length(A_range), cost_ratio_range(i), p.A, p.B);
    fprintf('=====================================================\n\n');
    
    % --- HYBRID OPTIMIZATION (QPSO + FMINCON) ---
    % 2.1 QPSO Global Search Phase
    fprintf('Starting QPSO global search phase...\n');
    num_control_vars = 2 * (N + 1);
    lb_controls = [zeros(1, N + 1), zeros(1, N + 1)];
    ub_controls = [p.umax * ones(1, N + 1), p.vmax * ones(1, N + 1)];
    particles = lb_controls + (ub_controls - lb_controls) .* rand(qpso.swarmSize, num_control_vars);
    pBest = particles; pBestCost = inf(qpso.swarmSize, 1);
    gBest = particles(1, :); gBestCost = inf;
    for iter = 1:qpso.maxIter
        parfor j = 1:qpso.swarmSize
            u = particles(j, 1:N+1); v = particles(j, N+2:end);
            cost = qpso_cost_func(u, v, T_guess, p, N);
            if cost < pBestCost(j), pBestCost(j) = cost; pBest(j, :) = particles(j, :); end
        end
        [min_pBestCost, idx] = min(pBestCost);
        if min_pBestCost < gBestCost, gBestCost = min_pBestCost; gBest = pBest(idx, :); end
        fprintf('  QPSO Iter: %d/%d | Best Cost: %.2f\n', iter, qpso.maxIter, gBestCost);
    end
    
    % 2.2 fmincon Local Refinement Phase
    fprintf('\nStarting fmincon refinement phase...\n');
    x_guess = linspace(p.x0, p.xf, N + 1); y_guess = p.y0 * ones(1, N + 1);
    z_guess = p.z0 * ones(1, N + 1); u_guess = gBest(1:N+1); v_guess = gBest(N+2:end);
    z0 = [x_guess, y_guess, z_guess, u_guess, v_guess, T_guess]';
    lb = [zeros(1, 3*(N+1)), lb_controls, 0.1]';
    ub = [inf(1, 3*(N+1)), ub_controls, inf]';
    options = optimoptions('fmincon', 'Algorithm', 'interior-point', 'Display', 'off', 'MaxIterations', 2000, 'MaxFunctionEvaluations', 5e5, 'ConstraintTolerance', 1e-6, 'OptimalityTolerance', 1e-6);
    [z_opt, J_opt] = fmincon(@(z) objective_func(z, p, N), z0, [], [], [], [], lb, ub, @(z) nonlinear_constraints(z, p, N), options);

    % --- Store Results for this iteration ---
    T_opt = z_opt(end);
    time_opt = linspace(0, T_opt, N+1);
    u_opt = z_opt((3*(N+1)+1):(4*(N+1)));
    v_opt = z_opt((4*(N+1)+1):(5*(N+1)));
    
    results.T_opt(i) = T_opt;
    results.J_opt(i) = J_opt;
    results.total_chemo(i) = trapz(time_opt, u_opt);
    results.total_immuno(i) = trapz(time_opt, v_opt);
    
end

%% 3. Post-Processing and Visualization of Sensitivity
fprintf('\n\n--- SENSITIVITY ANALYSIS COMPLETE ---\n');
fprintf('Plotting final sensitivity results...\n');

fig = figure('Name', 'Sensitivity to Relative Drug Cost (A/B)', 'Position', [100, 100, 900, 500]);
set(fig, 'Visible', 'on');
sgtitle('Sensitivity Analysis: Impact of Relative Drug Cost (A/B)', 'FontSize', 16, 'FontWeight', 'bold');

% Plot Total Drug Doses
subplot(1, 2, 1);
hold on;
plot(cost_ratio_range, results.total_chemo, 'k-o', 'LineWidth', 2, 'MarkerFaceColor', 'k', 'DisplayName', 'Total Chemo (u)');
plot(cost_ratio_range, results.total_immuno, 'm--s', 'LineWidth', 2, 'MarkerFaceColor', 'm', 'DisplayName', 'Total Immuno (v)');
hold off;
title('Total Drug Usage');
xlabel('Relative Drug Cost (A/B)');
ylabel('Total Integrated Dose');
grid on;
box on;
legend('show', 'Location', 'best');
set(gca, 'XScale', 'log'); % Use a log scale for the x-axis if ratios are far apart

% Plot Final Time
subplot(1, 2, 2);
plot(cost_ratio_range, results.T_opt, 'b-d', 'LineWidth', 2, 'MarkerFaceColor', 'b');
title('Optimal Treatment Duration');
xlabel('Relative Drug Cost (A/B)');
ylabel('Final Time T (days)');
grid on;
box on;
set(gca, 'XScale', 'log');

drawnow;

%% ==========================================================================
%                       LOCAL FUNCTION DEFINITIONS
%==========================================================================
% --- Objective and Constraint Functions for FMINCON ---
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
    
    fx = @(x,y,z,u,v,p) p.xi*x.*(1-x/p.x_inf) - p.theta*x.*y - p.alpha*x.*(u./(1+u));
    fy = @(x,y,z,u,v,p) p.a*(1-p.b*x).*y.*z + p.gamma - p.delta*y - p.kappa*y.*(u./(1+u)) + p.nu*y.*(v./(1+v));
    fz = @(x,y,z,u,v,p) p.sigma*x + p.psi*x.*(u./(1+u)) - p.mu*z;
    ceq_x = x(2:end) - x(1:end-1) - (h/2) * (fx(x(1:end-1), y(1:end-1), z_vals(1:end-1), u(1:end-1), v(1:end-1), p) + fx(x(2:end), y(2:end), z_vals(2:end), u(2:end), v(2:end), p));
    ceq_y = y(2:end) - y(1:end-1) - (h/2) * (fy(x(1:end-1), y(1:end-1), z_vals(1:end-1), u(1:end-1), v(1:end-1), p) + fy(x(2:end), y(2:end), z_vals(2:end), u(2:end), v(2:end), p));
    ceq_z = z_vals(2:end) - z_vals(1:end-1) - (h/2) * (fz(x(1:end-1), y(1:end-1), z_vals(1:end-1), u(1:end-1), v(1:end-1), p) + fz(x(2:end), y(2:end), z_vals(2:end), u(2:end), v(2:end), p));
    
    ceq_boundary = [x(1) - p.x0; y(1) - p.y0; z_vals(1) - p.z0; x(end) - p.xf];
    ceq = [ceq_x; ceq_y; ceq_z; ceq_boundary];
    c = [];
end
% --- Cost Function for QPSO (simpler, fixed T) ---
function cost = qpso_cost_func(u, v, T, p, N)
    time = linspace(0, T, N + 1);
    
    % Integrate the state dynamics
    w0 = [p.x0; p.y0; p.z0];
    try
        [~, w_traj] = ode45(@(t,w) dynamics(t,w,u,v,p,time), time, w0);
    catch
        cost = inf; % Return infinite cost if integration fails
        return;
    end
    
    x = w_traj(:,1)'; y = w_traj(:,2)'; z = w_traj(:,3)';
    
    % Check for NaN or Inf in states which indicates instability
    if any(isinf(x)) || any(isnan(x))
        cost = inf;
        return;
    end
    
    % Interpolate controls for integration
    u_interp = interp1(time, u, time, 'linear');
    v_interp = interp1(time, v, time, 'linear');
    
    % Calculate cost
    integrand = p.l1*x + p.l2*y + p.l3*z + p.A*u_interp + p.B*v_interp + p.C;
    cost = trapz(time, integrand);
    
    % Add a large penalty if the terminal constraint is not met
    penalty_weight = 1e6;
    cost = cost + penalty_weight * (x(end) - p.xf)^2;
end
function dwdt = dynamics(t, w, u_vec, v_vec, p, time_vec)
    % Unpack states
    x = w(1); y = w(2); z = w(3);
    
    % Interpolate controls to get value at current time t
    u = interp1(time_vec, u_vec, t, 'linear', 'extrap');
    v = interp1(time_vec, v_vec, t, 'linear', 'extrap');
    
    % System Dynamics
    dxdt = p.xi*x*(1-x/p.x_inf) - p.theta*x*y - p.alpha*x*(u/(1+u));
    dydt = p.a*(1-p.b*x)*y*z + p.gamma - p.delta*y - p.kappa*y*(u/(1+u)) + p.nu*y*(v/(1+v));
    dzdt = p.sigma*x + p.psi*x*(u/(1+u)) - p.mu*z;
    
    dwdt = [dxdt; dydt; dzdt];
end