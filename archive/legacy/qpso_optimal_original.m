%==========================================================================
% MATLAB Code for Solving and Verifying the Chemo-Immunotherapy OCP
% Method: Hybrid QPSO + fmincon (Direct Method)
% Model: Michaelis-Menten Kinetics for Controls
% Author: Bereket Sitotaw Kidane
% Status: Smoothing REMOVED to preserve sharp Bang-start
%==========================================================================

%% 1. Initialization and Setup
clear; clc; close all;
rng(1); % for reproducibility

fprintf('1. Setting up the optimal control problem with Michaelis-Menten kinetics...\n');

% --- Discretization and Time ---
N = 100; % Number of grid points for discretization
T_guess = 4.0; % Initial guess for final time T

% --- Model Parameters (Structure 'p') ---
p.xi = 0.5618;    p.x_inf = 780;    p.theta = 1.0;    p.alpha = 1.0;
p.a = 0.00726;    p.b = 0.00264;    p.gamma = 0.1181; p.delta = 0.3745;
p.kappa = 2.0;    p.nu = 1.0;       p.sigma = 0.1;    p.psi = 2.0;
p.mu = 0.45;

% --- Cost Function Weights (Scenario 1: Costly Chemo, Free Immuno) ---
p.l1 = 0;     p.l2 = 0;     p.l3 = 0; 
p.A = 1;      p.B = 0.3;    p.C = 1.0;

% --- Boundary Conditions and Control Limits ---
p.x0 = 650;      p.y0 = 0.1;      p.z0 = 400;      p.xf = 1.0;
p.umax = 1.0;    % Scaled Max chemotherapy dose (Saturation Test)
p.vmax = 1.0;    % Scaled Max immunotherapy dose (Saturation Test)

% --- QPSO Parameters ---
qpso.swarmSize = 150;
qpso.maxIter = 200;

%% 2. Run Hybrid Optimization (QPSO + fmincon)
fprintf('2. Starting Hybrid Optimization...\n');

% --- QPSO Global Search Phase ---
fprintf('   -> Phase 1: QPSO global search (fixed T = %.1f)...\n', T_guess);
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
    
    fprintf('      QPSO Iteration: %d / %d | Best Cost Found: %.4f\n', iter, qpso.maxIter, gBestCost);
    
    mBest = mean(pBest, 1);
    for i = 1:qpso.swarmSize
        phi = rand(1, num_control_vars);
        p_i = phi .* pBest(i, :) + (1 - phi) .* gBest;
        u_val = rand(1, num_control_vars);
        L = 2 * abs(particles(i, :) - mBest); 
        
        noise_amplitude = 0.1; 
        eta = randn(1, num_control_vars); 
        uniform_noise = rand(1, num_control_vars) - 0.5; 
        exploration_term = noise_amplitude * eta .* uniform_noise;
        
        particles(i, :) = p_i + (-1)^ceil(0.5+rand) .* L .* log(1 ./ u_val) + exploration_term;
    end
    particles = max(min(particles, ub_controls), lb_controls);
end
qpso_time = toc;
fprintf('   -> QPSO phase complete in %.2f seconds.\n\n', qpso_time);

% --- Direct Pass-Through (NO SMOOTHING) ---
fprintf('   -> Passing raw QPSO result to SQP (Preserving sharp transitions)...\n');
u_guess_raw = gBest(1:N+1);
v_guess_raw = gBest(N+2:end);

% Assign directly without smoothdata
u_guess_final = u_guess_raw;
v_guess_final = v_guess_raw;

% OPTIONAL SAFETY: If the start is high, clamp it to max to ensure clean "Bang" start
if u_guess_final(1) > 0.8 * p.umax
    u_guess_final(1) = p.umax;
end

% --- fmincon Local Refinement Phase ---
fprintf('   -> Phase 2: fmincon local refinement...\n');
tic;

x_guess = linspace(p.x0, p.xf, N + 1);
y_guess = p.y0 * ones(1, N + 1);
z_guess = p.z0 * ones(1, N + 1);

z0 = [x_guess, y_guess, z_guess, u_guess_final, v_guess_final, T_guess]';

lb_x = zeros(1, N + 1); ub_x = p.x_inf * ones(1,N+1);
lb_y = zeros(1, N + 1); ub_y = inf(1, N + 1);
lb_z = zeros(1, N + 1); ub_z = inf(1, N + 1);
lb_u = zeros(1, N + 1); ub_u = p.umax * ones(1, N + 1);
lb_v = zeros(1, N + 1); ub_v = p.vmax * ones(1, N + 1);
lb_T = 0.1;             ub_T = 20;

lb = [lb_x, lb_y, lb_z, lb_u, lb_v, lb_T]';
ub = [ub_x, ub_y, ub_z, ub_u, ub_v, ub_T]';

% --- Using SQP algorithm ---
options = optimoptions('fmincon', 'Algorithm', 'sqp', ...
    'Display', 'iter', 'MaxIterations', 5000, 'MaxFunctionEvaluations', 1e6, ...
    'ConstraintTolerance', 1e-6, 'OptimalityTolerance', 1e-6);

[z_opt, J_opt, ~, ~, lambda_struct]  = fmincon(@(z) objective_func(z, p, N), z0, ...
    [], [], [], [], lb, ub, @(z) nonlinear_constraints(z, p, N), options);

fmincon_time = toc;
fprintf('   -> fmincon refinement complete in %.2f seconds. Final cost J = %.4f\n\n', fmincon_time, J_opt);

%% 3. Post-Processing and Visualization of Optimal Solution
fprintf('3. Extracting and plotting the final optimal solution...\n');
idx = 1:(N+1);
x_opt = z_opt(idx);
y_opt = z_opt(idx + (N+1));
z_opt_vals = z_opt(idx + 2*(N+1));
u_opt = z_opt(idx + 3*(N+1));
v_opt = z_opt(idx + 4*(N+1));
T_opt = z_opt(end);
time_opt = linspace(0, T_opt, N + 1);

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');
set(groot, 'defaultAxesFontSize', 14);
plotLineWidth = 2.0;

figure('Name', 'Optimal Chemo-Immunotherapy Solution (Michaelis-Menten)', 'Position', [100, 100, 1200, 800]);
sgtitle(sprintf('Optimal Therapy Strategy (Final Time T = %.2f, Cost J = %.2f)', T_opt, J_opt), 'FontSize', 16, 'FontWeight', 'bold');

subplot(2, 2, 1); plot(time_opt, x_opt, 'r-', 'LineWidth', plotLineWidth); title('Tumor Volume ($x$)'); xlabel('Time'); ylabel('Population'); grid on; xlim([0, T_opt]);
subplot(2, 2, 2); plot(time_opt, y_opt, 'b-', 'LineWidth', plotLineWidth); title('Immune Cells ($y$)'); xlabel('Time'); ylabel('Density'); grid on; xlim([0, T_opt]);
subplot(2, 2, 3); plot(time_opt, z_opt_vals, 'g-', 'LineWidth', plotLineWidth); title('Tumor Antigen ($z$)'); xlabel('Time'); ylabel('Concentration'); grid on; xlim([0, T_opt]);

subplot(2, 2, 4); hold on; 
plot(time_opt, u_opt, 'k-', 'LineWidth', plotLineWidth, 'DisplayName', 'Chemo ($u$)'); 
plot(time_opt, v_opt, 'm--', 'LineWidth', plotLineWidth, 'DisplayName', 'Immuno ($v$)'); 
hold off; title('Optimal Control Dosages'); xlabel('Time'); ylabel('Dose'); grid on; legend('show', 'Location', 'best'); 
ylim([-0.1, max(p.umax, p.vmax) + 0.1]); 
xlim([0, T_opt]);

%% 4. PMP Optimality Verification
fprintf('4. Running PMP optimality verification on the solution...\n');

% --- Get Costates (Lagrange Multipliers) ---
lambda_multipliers = lambda_struct.eqnonlin;
h = T_opt / N; 
lambda1 = -lambda_multipliers(1:N) / h;
lambda2 = -lambda_multipliers(N+1 : 2*N) / h;
lambda3 = -lambda_multipliers(2*N+1 : 3*N) / h;

% Append last value to match time vector length
lambda1_opt = [lambda1; lambda1(end)];
lambda2_opt = [lambda2; lambda2(end)];
lambda3_opt = [lambda3; lambda3(end)];

% --- Calculate Switching Functions and Partials ---
switching_term_u = (-lambda1_opt .* p.alpha .* x_opt ...
                   - lambda2_opt .* p.kappa .* y_opt ...
                   + lambda3_opt .* p.psi   .* x_opt);
switching_term_v = lambda2_opt .* p.nu .* y_opt;

dH_du = p.A + switching_term_u ./ (1 + u_opt).^2;
dH_dv = p.B + switching_term_v ./ (1 + v_opt).^2;

d2H_du2 = -2 * switching_term_u ./ (1 + u_opt).^3;
d2H_dv2 = -2 * switching_term_v ./ (1 + v_opt).^3;

dt = T_opt / N;
dH_du_dt = gradient(dH_du, dt);
dH_dv_dt = gradient(dH_dv, dt);

% --- Calculate Hamiltonian Components ---
% 1. Control-dependent part of Hamiltonian
H_u_contribution = p.A * u_opt + switching_term_u .* (u_opt ./ (1 + u_opt));
H_v_contribution = p.B * v_opt + switching_term_v .* (v_opt ./ (1 + v_opt));
H_control_contribution = H_u_contribution + H_v_contribution;

% 2. State-dependent (drift) part of Hamiltonian
H_state_lagrangian = p.l1*x_opt + p.l2*y_opt + p.l3*z_opt_vals + p.C;
H_state_dynamics = lambda1_opt .* (p.xi*x_opt.*(1-x_opt/p.x_inf) - p.theta*x_opt.*y_opt) ...
                 + lambda2_opt .* (p.a*(1-p.b*x_opt).*y_opt.*z_opt_vals + p.gamma - p.delta*y_opt) ...
                 + lambda3_opt .* (p.sigma*x_opt - p.mu*z_opt_vals);
H_state_contribution = H_state_lagrangian + H_state_dynamics;

% 3. Total Hamiltonian (should be zero)
H_total = H_state_contribution + H_control_contribution;

% --- Generate PMP Verification Plots ---
figure('Name', 'PMP Optimality Verification (Michaelis-Menten)', 'Position', [200, 200, 1200, 800]);
sgtitle('PMP Optimality Conditions Analysis', 'FontSize', 16, 'FontWeight', 'bold');

% Subplot 1: Hamiltonian Components
subplot(2, 2, 1); hold on;
plot(time_opt, H_state_contribution, 'r:', 'LineWidth', 2, 'DisplayName', 'State Contribution');
plot(time_opt, H_control_contribution, 'b-', 'LineWidth', 2, 'DisplayName', 'Control Contribution');
plot(time_opt, H_total, 'k-', 'LineWidth', 3, 'DisplayName', 'Total Hamiltonian $H(t)$');
hold off; title('Hamiltonian Components (Free Time $\Rightarrow H \equiv 0$)'); 
xlabel('Time (years)'); ylabel('Value'); grid on; legend('show', 'Location', 'best'); xlim([0, T_opt]);
ylim(ylim + [-10, 10]); % Add padding to see the zero line

% Subplot 2: First Partials
subplot(2, 2, 2); hold on;
plot(time_opt, dH_du, 'b-', 'LineWidth', 2, 'DisplayName', '$\partial H / \partial u$');
plot(time_opt, dH_dv, 'c--', 'LineWidth', 2, 'DisplayName', '$\partial H / \partial v$');
plot(time_opt, zeros(size(time_opt)), 'k:', 'LineWidth', 1, 'HandleVisibility','off');
hold off; title('First Partials (Optimality Condition)'); xlabel('Time (years)'); ylabel('Value'); grid on; legend('show', 'Location', 'best'); xlim([0, T_opt]);

% Subplot 3: Second Partials
subplot(2, 2, 3); hold on;
plot(time_opt, d2H_du2, 'g-', 'LineWidth', 2, 'DisplayName', '$\partial^2H / \partial u^2$');
plot(time_opt, d2H_dv2, 'color', [0.3 0.9 0.3], 'LineStyle', '--', 'LineWidth', 2, 'DisplayName', '$\partial^2H / \partial v^2$');
plot(time_opt, zeros(size(time_opt)), 'k:', 'LineWidth', 1, 'HandleVisibility','off');
hold off; title('Second Partials (Legendre-Clebsch: $\partial^2H / \partial u^2 \ge 0$)'); 
xlabel('Time (years)'); ylabel('Value'); grid on; legend('show', 'Location', 'best'); xlim([0, T_opt]);

% Subplot 4: Time Derivatives of First Partials
subplot(2, 2, 4); hold on;
plot(time_opt, dH_du_dt, 'b-', 'LineWidth', 2, 'DisplayName', '$d/dt(\partial H / \partial u)$');
plot(time_opt, dH_dv_dt, 'c--', 'LineWidth', 2, 'DisplayName', '$d/dt(\partial H / \partial v)$');
plot(time_opt, zeros(size(time_opt)), 'k:', 'LineWidth', 1, 'HandleVisibility','off');
hold off; title('Time Derivatives of First Partials'); 
xlabel('Time (years)'); ylabel('Rate of Change'); grid on; legend('show', 'Location', 'best'); xlim([0, T_opt]);
%% 5. Performance Summary
fprintf('5. Generating Performance Summary Tables...\n\n');

% --- 1. Overall Performance ---
varNames = {'Metric', 'Value', 'Units'};
Metric = ["Total Therapy Duration"; "Final Objective Cost"];
Value = [T_opt; J_opt];
Units = ["years"; "N/A"];
T_Overall = table(Metric, Value, Units);

% --- 2. State Dynamics ---
[peak_y, idx_y] = max(y_opt);
[peak_z, idx_z] = max(z_opt_vals);
tumor_reduction_pct = (p.x0 - x_opt(end)) / p.x0 * 100;

Metric = ["Tumor Reduction"; "Initial Tumor"; "Final Tumor"; "Peak Immune Response (y)"; "Time of Peak (y)"; "Peak Antigen Level (z)"; "Time of Peak (z)"];
Value = [tumor_reduction_pct; p.x0; x_opt(end); peak_y; time_opt(idx_y); peak_z; time_opt(idx_z)];
Units = ["%"; "population"; "population"; "density"; "years"; "concentration"; "years"];
T_Dynamics = table(Metric, Value, Units);

% --- 3. Control Strategy Analysis ---
% Thresholds relative to dynamic Max Bounds (to correctly classify 1.5 as Bang)
BANG_OFF_THRESH = 0.01;

% Analyze Chemotherapy (u)
u_phase_bang_on = (u_opt > 0.99 * p.umax);
u_phase_singular = (u_opt >= BANG_OFF_THRESH & u_opt <= 0.99 * p.umax);
u_phase_bang_off = (u_opt < BANG_OFF_THRESH);

u_dur_on = trapz(time_opt, u_phase_bang_on);
u_dur_sing = trapz(time_opt, u_phase_singular);
u_dur_off = trapz(time_opt, u_phase_bang_off);
u_intensity = trapz(time_opt, u_opt) / T_opt * 100;

Metric = ["Bang-On Duration"; "Singular Arc Duration"; "Bang-Off Duration"; "Dosing Intensity (Avg)"];
Value = [u_dur_on; u_dur_sing; u_dur_off; u_intensity];
Units = ["years"; "years"; "years"; "%"];
T_Chemo = table(Metric, Value, Units);

% Analyze Immunotherapy (v)
v_phase_bang_on = (v_opt > 0.99 * p.vmax);
v_phase_singular = (v_opt >= BANG_OFF_THRESH & v_opt <= 0.99 * p.vmax);
v_phase_bang_off = (v_opt < BANG_OFF_THRESH);

v_dur_on = trapz(time_opt, v_phase_bang_on);
v_dur_sing = trapz(time_opt, v_phase_singular);
v_dur_off = trapz(time_opt, v_phase_bang_off);
v_intensity = trapz(time_opt, v_opt) / T_opt * 100;

Metric = ["Bang-On Duration"; "Singular Arc Duration"; "Bang-Off Duration"; "Dosing Intensity (Avg)"];
Value = [v_dur_on; v_dur_sing; v_dur_off; v_intensity];
Units = ["years"; "years"; "years"; "%"];
T_Immuno = table(Metric, Value, Units);

% --- 4. Computational Efficiency ---
Metric = ["QPSO Global Search"; "fmincon Local Refinement"; "Total Computation Time"];
Value = [qpso_time; fmincon_time; qpso_time + fmincon_time];
Units = ["seconds"; "seconds"; "seconds"];
T_Computation = table(Metric, Value, Units);

% --- Display Tables ---
fprintf('============================================================\n');
fprintf('           OPTIMAL THERAPY PERFORMANCE SUMMARY\n');
fprintf('============================================================\n\n');
fprintf('--- Overall Results ---\n');
disp(T_Overall);
fprintf('\n--- State Dynamics ---\n');
disp(T_Dynamics);
fprintf('\n--- Control Strategy: Chemotherapy (u) ---\n');
disp(T_Chemo);
fprintf('\n--- Control Strategy: Immunotherapy (v) ---\n');
disp(T_Immuno);
fprintf('\n--- Computational Efficiency ---\n');
disp(T_Computation);
fprintf('6. Script finished.\n');

%% ========================================================================
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
    [~, w_traj] = ode45(@(t,w) dynamics(t,w,u,v,p,time), time, w0, options);
    x = w_traj(:,1)'; y = w_traj(:,2)'; z_vals = w_traj(:,3)';
    
    integrand = p.l1*x + p.l2*y + p.l3*z_vals + p.A*u + p.B*v + p.C;
    cost = trapz(time, integrand);
    
    % CORRECTION: Penalize ANY deviation > 1.0
    cost = cost + 1e4 * (x(end) - p.xf)^2;
end

function dwdt = dynamics(t, w, u_vec, v_vec, p, time_vec)
    x = w(1); y = w(2); z = w(3);
    u = interp1(time_vec, u_vec, t, 'linear', 'extrap');
    v = interp1(time_vec, v_vec, t, 'linear', 'extrap');
    
    dxdt = p.xi*x*(1-x/p.x_inf) - p.theta*x*y - p.alpha*x*(u/(1+u));
    dydt = p.a*(1-p.b*x)*y*z + p.gamma - p.delta*y - p.kappa*y*(u/(1+u)) + p.nu*y*(v/(1+v));
    dzdt = p.sigma*x + p.psi*x*(u/(1+u)) - p.mu*z;
    
    dwdt = [dxdt; dydt; dzdt];
end