function plot_chemo_immuno_results(solution, summary, save_dir)
%PLOT_CHEMO_IMMUNO_RESULTS Paper-style state/control and PMP diagnostic figures.

    if nargin < 3
        save_dir = '';
    end

    set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
    set(groot, 'defaultLegendInterpreter', 'latex');
    set(groot, 'defaultTextInterpreter', 'latex');
    set(groot, 'defaultAxesFontSize', 14);
    lw = 2.0;

    t = solution.time;
    T_opt = solution.T;
    J_opt = solution.J;

    fig1 = figure('Name', 'Optimal Chemo-Immunotherapy Solution', 'Position', [100, 100, 1200, 800]);
    sgtitle(sprintf('%s (T = %.2f, J = %.2f)', solution.case_name, T_opt, J_opt), ...
        'FontSize', 16, 'FontWeight', 'bold');

    subplot(2, 2, 1);
    plot(t, solution.x, 'r-', 'LineWidth', lw);
    title('Tumor Volume ($x$)'); xlabel('Time (days)'); ylabel('Population'); grid on; xlim([0, T_opt]);

    subplot(2, 2, 2);
    plot(t, solution.y, 'b-', 'LineWidth', lw);
    title('Immune Cells ($y$)'); xlabel('Time (days)'); ylabel('Density'); grid on; xlim([0, T_opt]);

    subplot(2, 2, 3);
    plot(t, solution.z, 'g-', 'LineWidth', lw);
    title('Tumor Antigen ($z$)'); xlabel('Time (days)'); ylabel('Concentration'); grid on; xlim([0, T_opt]);

    subplot(2, 2, 4);
    hold on;
    plot(t, solution.u, 'k-', 'LineWidth', lw, 'DisplayName', 'Chemo ($u$)');
    plot(t, solution.v, 'm--', 'LineWidth', lw, 'DisplayName', 'Immuno ($v$)');
    hold off;
    title('Optimal Control Dosages'); xlabel('Time (days)'); ylabel('Dose'); grid on;
    legend('show', 'Location', 'best');
    ylim([-0.1, max(solution.p.umax, solution.p.vmax) + 0.1]);
    xlim([0, T_opt]);

    if ~isempty(solution.pmp)
        pmp = solution.pmp;
        fig2 = figure('Name', 'PMP Optimality Verification', 'Position', [200, 200, 1200, 800]);
        sgtitle('PMP Optimality Conditions Analysis', 'FontSize', 16, 'FontWeight', 'bold');

        subplot(2, 2, 1);
        hold on;
        plot(t, pmp.H_state_contribution, 'r:', 'LineWidth', 2, 'DisplayName', 'State Contribution');
        plot(t, pmp.H_control_contribution, 'b-', 'LineWidth', 2, 'DisplayName', 'Control Contribution');
        plot(t, pmp.H_total, 'k-', 'LineWidth', 3, 'DisplayName', 'Total Hamiltonian $H(t)$');
        hold off;
        title('Hamiltonian Components'); xlabel('Time (days)'); ylabel('Value');
        grid on; legend('show', 'Location', 'best'); xlim([0, T_opt]);

        subplot(2, 2, 2);
        hold on;
        plot(t, pmp.dH_du, 'b-', 'LineWidth', 2, 'DisplayName', '$\partial H / \partial u$');
        plot(t, pmp.dH_dv, 'c--', 'LineWidth', 2, 'DisplayName', '$\partial H / \partial v$');
        hold off;
        title('First Partials (Optimality Condition)'); xlabel('Time (days)'); ylabel('Value');
        grid on; legend('show', 'Location', 'best'); xlim([0, T_opt]);

        subplot(2, 2, 3);
        hold on;
        plot(t, pmp.d2H_du2, 'g-', 'LineWidth', 2, 'DisplayName', '$\partial^2H / \partial u^2$');
        plot(t, pmp.d2H_dv2, 'Color', [0.3 0.9 0.3], 'LineStyle', '--', 'LineWidth', 2, ...
            'DisplayName', '$\partial^2H / \partial v^2$');
        hold off;
        title('Second Partials (Legendre-Clebsch)'); xlabel('Time (days)'); ylabel('Value');
        grid on; legend('show', 'Location', 'best'); xlim([0, T_opt]);

        subplot(2, 2, 4);
        hold on;
        plot(t, pmp.dH_du_dt, 'b-', 'LineWidth', 2, 'DisplayName', '$d/dt(\partial H / \partial u)$');
        plot(t, pmp.dH_dv_dt, 'c--', 'LineWidth', 2, 'DisplayName', '$d/dt(\partial H / \partial v)$');
        hold off;
        title('Time Derivatives of First Partials'); xlabel('Time (days)'); ylabel('Rate of Change');
        grid on; legend('show', 'Location', 'best'); xlim([0, T_opt]);
    end

    if ~isempty(save_dir)
        if ~exist(save_dir, 'dir')
            mkdir(save_dir);
        end
        safe_name = matlab.lang.makeValidName(summary.case_name);
        exportgraphics(fig1, fullfile(save_dir, [safe_name '_states_controls.png']), 'Resolution', 200);
        if ~isempty(solution.pmp)
            exportgraphics(fig2, fullfile(save_dir, [safe_name '_pmp.png']), 'Resolution', 200);
        end
    end
end
