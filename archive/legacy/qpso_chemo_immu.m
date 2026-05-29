%==========================================================================
% QPSO-based Optimal Control for Chemo-Immunotherapy (Low-Dimensional)
%
% This script adapts the QPSO direct method to solve the optimal control
% problem for the 3-state model (x, y, z) from Ledzewicz et al. (2024).
% This direct method is more robust than the PMP shooting method.
%
%==========================================================================

%% 1. Configuration & Parameters
clear;
clc;
close all;
rng(42); % Set seed for reproducibility

fprintf('Setting up parameters for the Low-Dimensional Chemo-Immunotherapy problem...\n');

% --- Time and discretization ---
T = 3.5;                        % Final time (years)
n = 100;                        % Time divisions
time = linspace(0, T, n + 1);   % Time grid
dt = time(2) - time(1);         % Time step

% --- MODIFIED: Low-Dimensional Model Parameters (from Ledzewicz et al. 2024) ---
param = struct();
param.xi = 0.5618;    % Tumor growth parameter
param.x_inf = 780;    % Tumor carrying capacity
param.theta = 1;      % Immune-mediated tumor kill rate
param.alpha = 1;      % Chemo killing parameter on x
param.a = 0.00726;    % Immune stimulation by antigen
param.b = 0.00264;    % Tumor suppression of immune system
param.gamma_immune = 0.1181; % Constant influx of T-cells
param.delta = 0.3745; % Natural death rate of T-cells
param.kappa = 2;      % Chemo killing parameter on y
param.nu = 1;         % Immune boost efficacy
param.sigma = 0.1;    % Intrinsic immunogenicity
param.psi = 2;        % Therapy-induced immunogenicity
param.mu = 0.45;      % Antigen elimination

% --- MODIFIED: Cost Function Weights ---
param.A = 1.0;     % Weight for chemotherapy u
param.B = 0.3;     % Weight for immunotherapy v
param.C = 1.0;     % Weight for time
param.beta_x = 1.0;% Running cost on tumor burden

% --- Add time parameters to struct ---
param.time = time;
param.dt = dt;
param.n = n;

%% 2. Initial Conditions and Optimization Settings
% --- MODIFIED: Initial conditions for x, y, z ---
x0 = 650; % Initial tumor volume
y0 = 0.1; % Initial immunocompetent cells
z0 = 400; % Initial tumor antigen
initial_state = [x0; y0; z0];

% --- MODIFIED: Optimization for TWO controls (u and v) ---
num_vars = 2 * (n + 1); % We are optimizing for u_vec and v_vec concatenated
Umax = 1.0;
Vmax = 1.0;
lb = zeros(1, num_vars);
ub = [Umax * ones(1, n + 1), Vmax * ones(1, n + 1)];

% QPSO settings
qpso.swarmSize = 150;
qpso.maxIter = 600;

%% 3. QPSO Optimization Phase
fprintf('Starting QPSO optimization...\n');
particles = lb + (ub - lb) .* rand(qpso.swarmSize, num_vars);
pBest = particles;
pBestCost = inf(qpso.swarmSize, 1);
gBest = particles(1, :);
gBestCost = computeCost(gBest, param, initial_state);

% Optimization loop
for iter = 1:qpso.maxIter
    for i = 1:qpso.swarmSize
        currentTotalCost = computeCost(particles(i, :), param, initial_state);
        if currentTotalCost < pBestCost(i)
            pBestCost(i) = currentTotalCost;
            pBest(i, :) = particles(i, :);
        end
        if currentTotalCost < gBestCost
            gBestCost = currentTotalCost;
            gBest = particles(i, :);
        end
    end
    
    fprintf('QPSO Iteration: %d / %d | Best Cost Found: %.2f\n', iter, qpso.maxIter, gBestCost);
    
    p_mean = mean(pBest, 1);
    phi_rand = rand(qpso.swarmSize, num_vars);
    u_rand = rand(qpso.swarmSize, num_vars);
    u_rand(u_rand == 0) = eps;
    p = phi_rand .* pBest + (1 - phi_rand) .* gBest;
    L = (1./u_rand) .* log(1./u_rand);
    particles = p + L .* abs(particles - p_mean) .* randn(qpso.swarmSize, num_vars);
    particles = max(min(particles, ub), lb);
end
fprintf('QPSO phase complete. Best cost found: %.2f\n', gBestCost);

%% 4. fmincon Refinement Phase
fprintf('\nStarting fmincon refinement...\n');
options = optimoptions('fmincon', 'Algorithm', 'interior-point', 'Display', 'iter', 'MaxIterations', 1000, 'MaxFunctionEvaluations', 20000);
objectiveFun = @(uv_vec) computeCost(uv_vec, param, initial_state);
[optimal_controls_vec, finalCost] = fmincon(objectiveFun, gBest, [], [], [], [], lb, ub, [], options);
fprintf('fmincon refinement complete. Final cost: %.2f\n', finalCost);

%% 5. Post-Processing and Analysis
fprintf('\nComputing final trajectories and analysis...\n');
u_opt = optimal_controls_vec(1 : n+1);
v_opt = optimal_controls_vec(n+2 : end);
final_trajectory = eulerMethod(u_opt, v_opt, param, initial_state);
x_opt = final_trajectory(1, :);
y_opt = final_trajectory(2, :);
z_opt = final_trajectory(3, :);

%% 6. Visualization
fprintf('Generating plots...\n');
figure('Name', 'QPSO Chemo-Immunotherapy Results', 'Position', [100, 100, 1000, 800]);
sgtitle('QPSO Optimized Chemo-Immunotherapy', 'FontSize', 16, 'FontWeight', 'bold');

subplot(3, 2, 1);
plot(time, u_opt, 'r-', 'LineWidth', 2);
title('Optimal Control: Chemo u*(t)');
xlabel('Time (years)'); ylabel('Dosage'); grid on; ylim([-0.1, 1.1]);

subplot(3, 2, 2);
plot(time, v_opt, 'b-', 'LineWidth', 2);
title('Optimal Control: Immuno v*(t)');
xlabel('Time (years)'); ylabel('Dosage'); grid on; ylim([-0.1, 1.1]);

subplot(3, 2, 3);
plot(time, x_opt, 'k-', 'LineWidth', 2);
title('State: Tumor Volume x(t)');
xlabel('Time (years)'); ylabel('Volume'); grid on;

subplot(3, 2, 4);
plot(time, y_opt, 'g-', 'LineWidth', 2);
title('State: Immune Density y(t)');
xlabel('Time (years)'); ylabel('Density'); grid on;

subplot(3, 2, 5);
plot(time, z_opt, 'm-', 'LineWidth', 2);
title('State: Tumor Antigen z(t)');
xlabel('Time (years)'); ylabel('Level'); grid on;

%% 7. Local Functions

% --- MODIFIED: computeCost for the new problem ---
function totalCost = computeCost(uv_vec, param, initial_state)
    % Unpack the combined control vector
    u_vec = uv_vec(1 : param.n+1);
    v_vec = uv_vec(param.n+2 : end);
    
    % Simulate the system dynamics
    W = eulerMethod(u_vec, v_vec, param, initial_state);
    x_traj = W(1, :);
    
    % Calculate the total cost using the objective function
    running_cost = trapz(param.time, param.A*u_vec + param.B*v_vec + param.C + param.beta_x*x_traj);
    totalCost = running_cost; % No terminal cost in this formulation
end

% --- MODIFIED: dynamics function for the 3-state model ---
function dWdt = dynamics(W, u_val, v_val, param)
    x = W(1); y = W(2); z = W(3);
    
    dx = param.xi*x*(1 - x/param.x_inf) - param.theta*x*y - param.alpha*x*u_val;
    dy = param.a*(1-param.b*x)*y*z + param.gamma_immune - param.delta*y - param.kappa*y*u_val + param.nu*y*v_val;
    dz = param.sigma*x + param.psi*x*u_val - param.mu*z;
    
    dWdt = [dx; dy; dz];
end

% --- MODIFIED: eulerMethod for the 3-state, 2-control system ---
function W = eulerMethod(u_vec, v_vec, param, initial_state)
    W = zeros(3, length(param.time));
    W(:, 1) = initial_state;
    
    u_interp = @(t) interp1(param.time, u_vec, t, 'pchip');
    v_interp = @(t) interp1(param.time, v_vec, t, 'pchip');

    for k = 1:length(param.time) - 1
        t = param.time(k);
        dWdt = dynamics(W(:, k), u_interp(t), v_interp(t), param);
        W(:, k + 1) = W(:, k) + param.dt * dWdt;
    end
end