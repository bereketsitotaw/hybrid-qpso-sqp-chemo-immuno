function [solution, summary] = solve_chemo_immuno_ocp(cfg)
%SOLVE_CHEMO_IMMUNO_OCP Hybrid QPSO + SQP solver for chemo-immunotherapy OCP.
%
% Implements the two-stage pipeline described in Section 4 of the paper:
%   Stage 1: QPSO warm start with ode45 rollout and terminal penalty
%   Stage 2: Trapezoidal collocation + fmincon (SQP)
%
% Inputs:
%   cfg - struct with fields: p, N, T_guess, qpso, sqp, case_name (optional)
%
% Outputs:
%   solution - optimal trajectories and PMP diagnostics
%   summary    - timing, costs, and metadata

    p = cfg.p;
    N = cfg.N;
    T_guess = cfg.T_guess;
    qpso = cfg.qpso;
    sqp = cfg.sqp;

    if isfield(cfg, 'case_name')
        case_name = cfg.case_name;
    else
        case_name = 'unspecified';
    end

    fprintf('========================================\n');
    fprintf('Hybrid QPSO-SQP: %s\n', case_name);
    fprintf('========================================\n');

    fprintf('Stage 1: QPSO global search (T_guess = %.2f)...\n', T_guess);
    [gBest, gBestCost, qpso_time, qpso] = run_qpso_phase(p, N, T_guess, qpso);
    fprintf('   QPSO complete in %.2f s | best augmented cost = %.4f\n', qpso_time, gBestCost);

    fprintf('Stage 2: fmincon SQP collocation refinement...\n');
    [z_opt, J_opt, lambda_struct, fmincon_time, sqp] = run_sqp_refinement(p, N, gBest, T_guess, sqp);
    fprintf('   SQP complete in %.2f s | final cost J = %.4f\n', fmincon_time, J_opt);

    solution = unpack_solution(z_opt, N);
    solution.J = J_opt;
    solution.p = p;
    solution.case_name = case_name;

    if cfg.run_pmp
        solution.pmp = compute_pmp_diagnostics( ...
            solution.x, solution.y, solution.z, solution.u, solution.v, ...
            solution.T, N, p, lambda_struct);
    else
        solution.pmp = [];
    end

    summary = struct();
    summary.case_name = case_name;
    summary.qpso_time = qpso_time;
    summary.fmincon_time = fmincon_time;
    summary.total_time = qpso_time + fmincon_time;
    summary.qpso_best_cost = gBestCost;
    summary.final_cost = J_opt;
    summary.T_final = solution.T;
    summary.total_chemo = trapz(solution.time, solution.u);
    summary.total_immuno = trapz(solution.time, solution.v);
    summary.terminal_tumor = solution.x(end);
    summary.cfg = cfg;

    fprintf('----------------------------------------\n');
    fprintf('T* = %.4f | J* = %.4f | x(T) = %.4f\n', summary.T_final, summary.final_cost, summary.terminal_tumor);
    fprintf('Total chemo = %.4f | Total immuno = %.4f\n', summary.total_chemo, summary.total_immuno);
    fprintf('========================================\n');
end
