%==========================================================================
% MATLAB Code for Solving and Verifying the Chemo-Immunotherapy OCP
% Method: Hybrid QPSO + fmincon (Direct Method)
% Model: Control-Affine (Bang-Bang)
% Author: Bereket Sitotaw Kidane
%==========================================================================
%% 1. Initialization and Setup
clear; clc; close all;
rng(1); % for reproducibility
fprintf('1. Setting up the optimal control problem...\n');

% --- Discretization and Time ---
N = 100; % Number of grid points for discretization
T_guess = 3.0; % Fixed time horizon for QPSO phase, fmincon will optimize it

% --- Model Parameters (from Ledzewicz et al., 2024, Table 1) ---
p.xi = 0.5618;    % Tumor growth rate
p.x_inf = 780;    % Tumor carrying capacity
p.theta = 1.0;    % Immune-mediated tumor elimination rate
p.alpha = 1.0;    % Chemotherapy efficacy on tumor
p.a = 0.00726;    % Immune stimulation rate by antigen
p.b = 0.00264;    % Tumor-induced immune suppression threshold
p.gamma = 0.1181; % Influx of immune cells
p.delta = 0.3745; % Natural death rate of immune cells
p.kappa = 2.0;    % Chemotherapy toxicity on immune cells
p.nu = 1.0;       % Immunotherapy efficacy
p.sigma = 0.1;    % Intrinsic antigen production rate
p.psi = 2.0;      % Therapy-induced antigen production rate
p.mu = 0.45;      % Antigen elimination rate

% --- Cost Function Weights (from Ledzewicz et al., 2024, Sec 5.1) ---
p.l1 = 0; p.l2 = 0; p.l3 = 0; % No running cost on states
p.A = 0.0;     % Cost of chemotherapy
p.B = 0.0;     % Cost of immunotherapy
p.C = 1.0;     % Cost of time (to minimize T)

% --- Boundary Conditions and Control Limits ---
p.x0 = 650;      % Initial tumor volume
p.y0 = 0.1;      % Initial immune cell density
p.z0 = 400;      % Initial antigen concentration
p.xf = 1.0;      % Target final tumor volume
p.umax = 3.0;    % Max chemotherapy dose
p.vmax = 3.0;    % Max immunotherapy dose

% --- QPSO Parameters ---
qpso.swarmSize = 150;
qpso.maxIter = 150; % Number of iterations for global search

%% 2. QPSO Global Search Phase
fprintf('2. Starting QPSO global search phase (fixed T = %.1f)...\n', T_guess);
tic;
num_control_vars = 2 * (N + 1);
lb_controls = [zeros(1, N + 1), zeros(1, N + 1)];
ub_controls = [p.umax * ones(1, N + 1), p.vmax * ones(1, N + 1)];
particles = lb_controls + (ub_controls - lb_controls) .* rand(qpso.swarmSize, num_control_vars);
pBest = particles;
pBestCost = inf(qpso.swarmSize, 1);
gBest = particles(1, :);
gBestCost = inf;

for iter = 1:qpso.maxIter
    parfor i = 1:qpso.swarmSize
        u = particles(i, 1:N+1);
        v = particles(i, N+2:end);
        cost = qpso_cost_func(u, v, T_guess, p, N);
        if cost < pBestCost(i)
            pBestCost(i) = cost;
            pBest(i, :) = particles(i, :);
        end
    end
    [min_pBestCost, idx] = min(pBestCost);
    if min_pBestCost < gBestCost
        gBestCost = min_pBestCost;
        gBest = pBest(idx, :);
    end
    fprintf('   QPSO Iteration: %d / %d | Best Cost Found: %.4f\n', iter, qpso.maxIter, gBestCost);
    mBest = mean(pBest, 1);
    for i = 1:qpso.swarmSize
        phi = rand(1, num_control_vars);
        p_i = phi .* pBest(i, :) + (1 - phi) .* gBest;
        u_val = rand(1, num_control_vars);
        L = (1 / norm(u_val)) * abs(particles(i, :) - p_i);
        particles(i, :) = p_i + (-1)^ceil(0.5+rand) .* L .* log(1 ./ u_val);
    end
    particles = max(min(particles, ub_controls), lb_controls);
end
qpso_time = toc;
fprintf('   QPSO phase complete in %.2f seconds.\n\n', qpso_time);

%% 3. fmincon Local Refinement Phase
fprintf('3. Starting fmincon refinement phase...\n');
tic;
x_guess = linspace(p.x0, p.xf, N + 1);
y_guess = p.y0 * ones(1, N + 1);
z_guess = p.z0 * ones(1, N + 1);
u_guess = gBest(1:N+1);
v_guess = gBest(N+2:end);
z0 = [x_guess, y_guess, z_guess, u_guess, v_guess, T_guess]';

lb_x = zeros(1, N + 1); ub_x = inf(1, N + 1);
lb_y = zeros(1, N + 1); ub_y = inf(1, N + 1);
lb_z = zeros(1, N + 1); ub_z = inf(1, N + 1);
lb_u = zeros(1, N + 1); ub_u = p.umax * ones(1, N + 1);
lb_v = zeros(1, N + 1); ub_v = p.vmax * ones(1, N + 1);
lb_T = 0.1;             ub_T = inf;
lb = [lb_x, lb_y, lb_z, lb_u, lb_v, lb_T]';
ub = [ub_x, ub_y, ub_z, ub_u, ub_v, ub_T]';

options = optimoptions('fmincon', 'Algorithm', 'interior-point', ...
    'Display', 'iter', 'MaxIterations', 5000, ...
    'MaxFunctionEvaluations', 1e6, 'ConstraintTolerance', 1e-6, ...
    'OptimalityTolerance', 1e-6);

[z_opt, J_opt, ~, ~, lambda_struct] = fmincon(@(z) objective_func(z, p, N), z0, ...
    [], [], [], [], lb, ub, @(z) nonlinear_constraints(z, p, N), options);
fmincon_time = toc;
fprintf('   fmincon refinement complete in %.2f seconds. Final cost J = %.4f\n\n', fmincon_time, J_opt);

%% 4. Post-Processing and Visualization
fprintf('4. Extracting and plotting the final optimal solution...\n');
N_pts = N + 1;
idx = 1:N_pts;
x_opt = z_opt(idx);
y_opt = z_opt(idx + N_pts);
z_opt_vals = z_opt(idx + 2*N_pts);
u_opt = z_opt(idx + 3*N_pts);
v_opt = z_opt(idx + 4*N_pts);
T_opt = z_opt(end);
time_opt = linspace(0, T_opt, N_pts);

figure('Name', 'Optimal Chemo-Immunotherapy Solution', 'Position', [100, 100, 1000, 600]);
sgtitle(sprintf('Optimal Therapy Strategy (Final Time T = %.3f)', T_opt), 'FontSize', 16, 'FontWeight', 'bold');
subplot(2, 2, 1); plot(time_opt, x_opt, 'r-', 'LineWidth', 2); title('Tumor Volume (x)'); xlabel('Time (days)'); ylabel('Population'); grid on; xlim([0, T_opt]);
subplot(2, 2, 2); plot(time_opt, y_opt, 'b-', 'LineWidth', 2); title('Immune Cells (y)'); xlabel('Time (days)'); ylabel('Density'); grid on; xlim([0, T_opt]);
subplot(2, 2, 3); plot(time_opt, z_opt_vals, 'g-', 'LineWidth', 2); title('Tumor Antigen (z)'); xlabel('Time (days)'); ylabel('Concentration'); grid on; xlim([0, T_opt]);
subplot(2, 2, 4); hold on;
plot(time_opt, u_opt, 'k-', 'LineWidth', 2, 'DisplayName', 'Chemo (u)');
plot(time_opt, v_opt, 'm--', 'LineWidth', 2, 'DisplayName', 'Immuno (v)');
hold off; title('Optimal Control Dosages'); xlabel('Time (days)'); ylabel('Dose'); grid on; legend('show', 'Location', 'best'); xlim([0, T_opt]);
% --- MODIFIED LINE ---
ylim([-0.1, max(p.umax, p.vmax) + 0.1]);

%% 5. PMP Optimality Verification
fprintf('5. Running PMP optimality verification on the solution...\n');
% --- Step 1: Recover Costates ---
lambda_multipliers = lambda_struct.eqnonlin;
h = T_opt / N; 
lambda1 = -lambda_multipliers(1:N) / h;
lambda2 = -lambda_multipliers(N+1 : 2*N) / h;
lambda3 = -lambda_multipliers(2*N+1 : 3*N) / h;
lambda1_opt = [lambda1; lambda1(end)];
lambda2_opt = [lambda2; lambda2(end)];
lambda3_opt = [lambda3; lambda3(end)];

% --- Step 2: Calculate All Necessary Hamiltonian Components ---
Phi_1 = p.A - lambda1_opt .* p.alpha .* x_opt ...
        - lambda2_opt .* p.kappa .* y_opt ...
        + lambda3_opt .* p.psi   .* x_opt;
Phi_2 = p.B + lambda2_opt .* p.nu .* y_opt;
Hu_t = Phi_1 .* u_opt;
Hv_t = Phi_2 .* v_opt;
d2H_du2 = zeros(size(time_opt));
d2H_dv2 = zeros(size(time_opt));
dt = T_opt / N;
dPhi1_dt = gradient(Phi_1, dt);
dPhi2_dt = gradient(Phi_2, dt);

% --- Step 3: Plot the 2x2 Optimality Dashboard ---
figure('Name', 'PMP Optimality Verification', 'Position', [200, 200, 1200, 800]);
sgtitle('PMP Optimality Conditions Analysis', 'FontSize', 16, 'FontWeight', 'bold');

subplot(2, 2, 1); hold on;
plot(time_opt, Hu_t, 'b-', 'LineWidth', 2, 'DisplayName', 'Chemo Term ($\Phi_1 u^*$)');
plot(time_opt, Hv_t, 'c--', 'LineWidth', 2, 'DisplayName', 'Immuno Term ($\Phi_2 v^*$)');
hold off;
title('Control Contribution to Hamiltonian'); xlabel('Time (days)'); ylabel('Value'); grid on; legend('show', 'Location', 'best'); xlim([0, T_opt]);

subplot(2, 2, 2); hold on;
plot(time_opt, Phi_1, 'b-', 'LineWidth', 2, 'DisplayName', '$\Phi_1$ for Chemo');
plot(time_opt, Phi_2, 'c--', 'LineWidth', 2, 'DisplayName', '$\Phi_2$ for Immuno');
plot(time_opt, zeros(size(time_opt)), 'k:', 'LineWidth', 1, 'HandleVisibility','off');
hold off;
title('First Partials (Switching Functions)'); xlabel('Time (days)'); ylabel('Value'); grid on; legend('show', 'Location', 'best'); xlim([0, T_opt]);

subplot(2, 2, 3); hold on;
plot(time_opt, d2H_du2, 'g-', 'LineWidth', 2, 'DisplayName', '$\partial^2H / \partial u^2$');
plot(time_opt, d2H_dv2, 'color', [0.3 0.9 0.3], 'LineStyle', '--', 'LineWidth', 2, 'DisplayName', '$\partial^2H / \partial v^2$');
hold off;
title('Second Partials (Zero for Control-Affine System)'); xlabel('Time (days)'); ylabel('Value'); grid on; legend('show', 'Location', 'best'); xlim([0, T_opt]); ylim([-1, 1]);

subplot(2, 2, 4); hold on;
plot(time_opt, dPhi1_dt, 'b-', 'LineWidth', 2, 'DisplayName', '$d\Phi_1/dt$');
plot(time_opt, dPhi2_dt, 'c--', 'LineWidth', 2, 'DisplayName', '$d\Phi_2/dt$');
plot(time_opt, zeros(size(time_opt)), 'k:', 'LineWidth', 1, 'HandleVisibility','off');
hold off;
title('Time Derivatives of Switching Functions'); xlabel('Time (days)'); ylabel('Rate of Change'); grid on; legend('show', 'Location', 'best'); xlim([0, T_opt]);

fprintf('6. Script finished.\n');

%% ==========================================================================
%                       LOCAL FUNCTION DEFINITIONS
%==========================================================================
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
    
    fx = @(x,y,z,u,v,p) p.xi*x.*(1-x/p.x_inf) - p.theta*x.*y - p.alpha*x.*u;
    fy = @(x,y,z,u,v,p) p.a*(1-p.b*x).*y.*z + p.gamma - p.delta*y - p.kappa*y.*u + p.nu*y.*v;
    fz = @(x,y,z,u,v,p) p.sigma*x + p.psi*x.*u - p.mu*z;
    
    ceq_x = x(2:end) - x(1:end-1) - (h/2) * (fx(x(1:end-1), y(1:end-1), z_vals(1:end-1), u(1:end-1), v(1:end-1), p) + fx(x(2:end), y(2:end), z_vals(2:end), u(2:end), v(2:end), p));
    ceq_y = y(2:end) - y(1:end-1) - (h/2) * (fy(x(1:end-1), y(1:end-1), z_vals(1:end-1), u(1:end-1), v(1:end-1), p) + fy(x(2:end), y(2:end), z_vals(2:end), u(2:end), v(2:end), p));
    ceq_z = z_vals(2:end) - z_vals(1:end-1) - (h/2) * (fz(x(1:end-1), y(1:end-1), z_vals(1:end-1), u(1:end-1), v(1:end-1), p) + fz(x(2:end), y(2:end), z_vals(2:end), u(2:end), v(2:end), p));
    
    ceq_boundary = [x(1) - p.x0; y(1) - p.y0; z_vals(1) - p.z0; x(end) - p.xf];
    
    ceq = [ceq_x; ceq_y; ceq_z; ceq_boundary];
    c = [];
end

function cost = qpso_cost_func(u, v, T, p, N)
    time = linspace(0, T, N + 1);
    w0 = [p.x0; p.y0; p.z0];
    [~, w_traj] = ode45(@(t,w) dynamics(t,w,u,v,p,time), time, w0);
    x = w_traj(:,1)'; y = w_traj(:,2)'; z = w_traj(:,3)';
    integrand = p.l1*x + p.l2*y + p.l3*z + p.A*u + p.B*v + p.C;
    cost = trapz(time, integrand);
    if abs(x(end) - p.xf) > 1.0 
        cost = cost + 1e6 * (x(end) - p.xf)^2;
    end
end

function dwdt = dynamics(t, w, u_vec, v_vec, p, time_vec)
    x = w(1); y = w(2); z = w(3);
    u = interp1(time_vec, u_vec, t, 'linear', 'extrap');
    v = interp1(time_vec, v_vec, t, 'linear', 'extrap');
    
    dxdt = p.xi*x*(1-x/p.x_inf) - p.theta*x*y - p.alpha*x*u;
    dydt = p.a*(1-p.b*x)*y*z + p.gamma - p.delta*y - p.kappa*y*u + p.nu*y*v;
    dzdt = p.sigma*x + p.psi*x*u - p.mu*z;
    
    dwdt = [dxdt; dydt; dzdt];
end