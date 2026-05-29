function plot_quadrotor_results(solution, summary, problem, opts)
%==========================================================================
% plot_quadrotor_results.m  -- FIXED DROP-IN
%==========================================================================
% Visualization utilities for the quadrotor OCP:
%   - Static trajectory figure (paper-ready)
%   - Controls figure
%   - Animation (optional MP4/GIF)
%   - PMP plots (optional)
%
% FIXES:
%   - Uses Hermite dense clearance check (x,y from vx,vy)
%   - Summary text supports BOTH old and new summary field names
%   - draw_circle uses proper EdgeColor + LineStyle (no 'k:' misuse)
%==========================================================================

    set(groot, 'defaultAxesTickLabelInterpreter','latex');
    set(groot, 'defaultLegendInterpreter','latex');
    set(groot, 'defaultTextInterpreter','latex');
    set(groot, 'defaultAxesFontSize', 14);

    if nargin < 4, opts = struct(); end
    opts = set_defaults(opts, struct( ...
        'show_reference', true, ...
        'show_arrows', true, ...
        'dense_check_sub', 12, ...
        'animate', true, ...
        'frame_skip', 1, ...
        'fps', 30, ...
        'save_mp4', false, ...
        'mp4_name', 'quad_traj.mp4', ...
        'save_gif', false, ...
        'gif_name', 'quad_traj.gif', ...
        'plot_pmp', true ...
    ));

    plot_trajectory_static(solution, summary, problem, opts);
    plot_controls(solution, problem);

    if opts.animate
        animate_trajectory(solution, problem, opts);
    end

    if opts.plot_pmp && isfield(solution,'pmp') && ~isempty(solution.pmp) && isfield(solution.pmp,'hamiltonian')
        plot_pmp(solution);
    end
end

%==========================================================================
% Static trajectory figure
%==========================================================================
function plot_trajectory_static(solution, summary, problem, opts)
    p = problem.params;
    t = solution.time_vec(:)';
    x = solution.states(1,:); y = solution.states(2,:);
    vx = solution.states(3,:); vy = solution.states(4,:);

    [xS,yS,xG,yG] = get_start_goal(solution, problem);

    [minClearDense, violMaskDense, xD, yD] = dense_clearance_check(solution, problem, opts.dense_check_sub);

    % Use map bounds if present
    if isfield(problem.bounds,'xlim') && isfield(problem.bounds,'ylim')
        xLim = problem.bounds.xlim; yLim = problem.bounds.ylim;
    elseif isfield(p,'map') && isfield(p.map,'xlim') && isfield(p.map,'ylim')
        xLim = p.map.xlim; yLim = p.map.ylim;
    else
        pad = 0.8;
        xLim = [min(x)-pad, max(x)+pad];
        yLim = [min(y)-pad, max(y)+pad];
    end

    figure('Name','Trajectory (Static)','Position',[80 80 940 680]);
    ax = axes; hold(ax,'on'); grid(ax,'on'); axis(ax,'equal');
    xlim(ax,xLim); ylim(ax,yLim);
    xlabel(ax,'X Position'); ylabel(ax,'Y Position');

    title(ax, sprintf('Trajectory | T=%.2f, J=%.2f | minClear(dense)=%.3f', ...
        solution.T_final, solution.cost, minClearDense), 'Interpreter','none');

    draw_obstacles(ax, p, t);

    if opts.show_reference
        plot(ax,[xS xG],[yS yG],'k:','LineWidth',1.2,'DisplayName','Straight-line ref');
    end

    plot(ax, x, y, 'b-','LineWidth',2.6,'DisplayName','Trajectory');

    if opts.show_arrows
        nArrow = 7;
        idx = unique(round(linspace(1, numel(x)-1, nArrow)));
        quiver(ax, x(idx), y(idx), x(idx+1)-x(idx), y(idx+1)-y(idx), 0, ...
            'LineWidth', 1.0, 'MaxHeadSize', 0.8, 'HandleVisibility','off');
    end

    plot(ax, xS, yS, 'g*','MarkerSize',12,'LineWidth',1.7,'DisplayName','Start');
    plot(ax, xG, yG, 'ro','MarkerSize',9,'LineWidth',1.7,'DisplayName','Goal');

    if any(violMaskDense)
        plot(ax, xD(violMaskDense), yD(violMaskDense), 'rx', 'MarkerSize',8,'LineWidth',1.5, ...
            'DisplayName','Dense collision points');
    end

    % ---- summary text: supports new and old field names ----
    [maxIneq, maxEq, constrV] = extract_solver_feas_text(summary);
    txt = sprintf('Solver feas: node max(max(0,c))=%.1e, node max|ceq|=%.1e, constr=%.1e', ...
        maxIneq, maxEq, constrV);

    text(ax, xLim(1)+0.02*range(xLim), yLim(2)-0.05*range(yLim), txt, ...
        'Interpreter','none','FontSize',11);

    % Optional: show max speed for sanity
    vnorm = sqrt(vx.^2 + vy.^2);
    if isfield(problem.bounds,'v_max')
        vmax = problem.bounds.v_max;
    else
        vmax = getfield_default(problem.params,'v_max',NaN);
    end
    txt2 = sprintf('Max ||v|| = %.2f (limit %.2f)', max(vnorm), vmax);
    text(ax, xLim(1)+0.02*range(xLim), yLim(2)-0.10*range(yLim), txt2, ...
        'Interpreter','none','FontSize',11);

    legend(ax,'show','Location','bestoutside');
    hold(ax,'off');
end

%==========================================================================
% Controls
%==========================================================================
function plot_controls(solution, problem)
    t = solution.time_vec(:)';
    U = solution.controls;

    figure('Name','Controls','Position',[120 120 900 420]);
    ax = axes; hold(ax,'on'); grid(ax,'on');
    title(ax,'Control Trajectories'); xlabel(ax,'Time'); ylabel(ax,'Acceleration');

    plot(ax, t, U(1,:), 'LineWidth',2.0,'DisplayName','$a_x$');
    plot(ax, t, U(2,:), 'LineWidth',2.0,'DisplayName','$a_y$');

    if isfield(problem.bounds,'control_min') && isfield(problem.bounds,'control_max')
        umin = problem.bounds.control_min(:);
        umax = problem.bounds.control_max(:);
        plot(ax, t, umin(1)*ones(size(t)), 'k:','HandleVisibility','off');
        plot(ax, t, umax(1)*ones(size(t)), 'k:','HandleVisibility','off');
        plot(ax, t, umin(2)*ones(size(t)), 'k--','HandleVisibility','off');
        plot(ax, t, umax(2)*ones(size(t)), 'k--','HandleVisibility','off');
        plot(ax, nan, nan, 'k:','DisplayName','$a_x$ bounds');
        plot(ax, nan, nan, 'k--','DisplayName','$a_y$ bounds');
    end

    legend(ax,'show','Location','best');
    hold(ax,'off');
end

%==========================================================================
% Animation
%==========================================================================
function animate_trajectory(solution, problem, opts)
    p = problem.params;
    t = solution.time_vec(:)';
    x = solution.states(1,:); y = solution.states(2,:);

    [xS,yS,xG,yG] = get_start_goal(solution, problem);

    if isfield(problem.bounds,'xlim') && isfield(problem.bounds,'ylim')
        xLim = problem.bounds.xlim; yLim = problem.bounds.ylim;
    elseif isfield(p,'map') && isfield(p.map,'xlim') && isfield(p.map,'ylim')
        xLim = p.map.xlim; yLim = p.map.ylim;
    else
        pad = 0.8;
        xLim = [min(x)-pad, max(x)+pad];
        yLim = [min(y)-pad, max(y)+pad];
    end

    fig = figure('Name','Trajectory Animation','Position',[140 140 940 680]);
    ax = axes(fig); hold(ax,'on'); grid(ax,'on'); axis(ax,'equal');
    xlim(ax,xLim); ylim(ax,yLim);
    xlabel(ax,'X Position'); ylabel(ax,'Y Position');

    draw_obstacles(ax, p, t);
    plot(ax, xS, yS, 'g*','MarkerSize',12,'LineWidth',1.7,'HandleVisibility','off');
    plot(ax, xG, yG, 'ro','MarkerSize',9,'LineWidth',1.7,'HandleVisibility','off');

    trajLine = plot(ax, nan, nan, 'b-','LineWidth',2.6);
    quadDot  = plot(ax, nan, nan, 'bo','MarkerSize',7,'LineWidth',1.5);

    % Writers
    if opts.save_mp4
        vw = VideoWriter(opts.mp4_name, 'MPEG-4');
        vw.FrameRate = opts.fps;
        open(vw);
    else
        vw = [];
    end
    if opts.save_gif
        if exist(opts.gif_name,'file'), delete(opts.gif_name); end
    end

    kList = 1:opts.frame_skip:numel(t);
    for kk = kList
        set(trajLine,'XData',x(1:kk),'YData',y(1:kk));
        set(quadDot, 'XData',x(kk),  'YData',y(kk));

        title(ax, sprintf('t=%.2f/%.2f | J=%.2f', t(kk), t(end), solution.cost), 'Interpreter','none');
        drawnow;

        if opts.save_mp4 || opts.save_gif
            fr = getframe(fig);
        end
        if opts.save_mp4
            writeVideo(vw, fr);
        end
        if opts.save_gif
            [im, cm] = rgb2ind(fr.cdata, 256);
            if kk == kList(1)
                imwrite(im, cm, opts.gif_name, 'gif', 'LoopCount', inf, 'DelayTime', 1/opts.fps);
            else
                imwrite(im, cm, opts.gif_name, 'gif', 'WriteMode', 'append', 'DelayTime', 1/opts.fps);
            end
        end
    end

    if opts.save_mp4, close(vw); fprintf('Saved MP4: %s\n', opts.mp4_name); end
    if opts.save_gif, fprintf('Saved GIF: %s\n', opts.gif_name); end

    hold(ax,'off');
end

%==========================================================================
% PMP plots
%==========================================================================
function plot_pmp(solution)
    t = solution.time_vec(:)';

    figure('Name','PMP Diagnostics','Position',[200 120 1100 720]);
    sgtitle('PMP / KKT Optimality Diagnostics','FontSize',16,'FontWeight','bold');

    subplot(2,2,1);
    plot(t, solution.pmp.hamiltonian,'k-','LineWidth',2); grid on;
    title('Hamiltonian trace'); xlabel('Time'); ylabel('H(t)');

    subplot(2,2,2);
    if isfield(solution.pmp,'dH_du') && ~isempty(solution.pmp.dH_du)
        plot(t, solution.pmp.dH_du(1,:), 'LineWidth',2,'DisplayName','$\partial H/\partial a_x$'); hold on;
        plot(t, solution.pmp.dH_du(2,:), 'LineWidth',2,'DisplayName','$\partial H/\partial a_y$'); grid on;
        yline(0,'k:'); title('Stationarity residual'); xlabel('Time'); ylabel('dH/du');
        legend('show'); hold off;
    else
        text(0.5,0.5,'dH/du not available','HorizontalAlignment','center'); axis off;
    end

    subplot(2,2,3);
    if isfield(solution.pmp,'d2H_du2') && ~isempty(solution.pmp.d2H_du2)
        plot(t, solution.pmp.d2H_du2(1,:), 'LineWidth',2,'DisplayName','$\partial^2H/\partial a_x^2$'); hold on;
        plot(t, solution.pmp.d2H_du2(2,:), 'LineWidth',2,'DisplayName','$\partial^2H/\partial a_y^2$'); grid on;
        title('Second partials'); xlabel('Time'); ylabel('d^2H/du^2');
        legend('show'); hold off;
    else
        text(0.5,0.5,'d2H/du2 not available','HorizontalAlignment','center'); axis off;
    end

    subplot(2,2,4);
    if isfield(solution.pmp,'dH_du_dt') && ~isempty(solution.pmp.dH_du_dt)
        plot(t, solution.pmp.dH_du_dt(1,:), 'LineWidth',2,'DisplayName','$d/dt(\partial H/\partial a_x)$'); hold on;
        plot(t, solution.pmp.dH_du_dt(2,:), 'LineWidth',2,'DisplayName','$d/dt(\partial H/\partial a_y)$'); grid on;
        title('Time derivative'); xlabel('Time'); ylabel('d/dt(dH/du)');
        legend('show'); hold off;
    else
        text(0.5,0.5,'d/dt(dH/du) not available','HorizontalAlignment','center'); axis off;
    end
end

%==========================================================================
% Drawing + dense collision check (Hermite)
%==========================================================================
function draw_obstacles(ax, p, t)
    if ~isfield(p,'obs') || isempty(p.obs), return; end

    for k = 1:numel(p.obs)
        ox = p.obs(k).pos0(1) + p.obs(k).vel(1)*t;
        oy = p.obs(k).pos0(2) + p.obs(k).vel(2)*t;
        Rk = p.obs(k).R + p.safety_margin;

        if norm(p.obs(k).vel) > 0
            plot(ax, ox, oy, 'k--','LineWidth',1.0,'HandleVisibility','off');
        end

        draw_circle(ax, ox(1),   oy(1),   Rk, ':',  1.4);   % keep-out at t=0
        draw_circle(ax, ox(end), oy(end), Rk, '-',  1.2);   % keep-out at t=T

        plot(ax, ox(1),   oy(1),   'kx','MarkerSize',6,'LineWidth',1.0,'HandleVisibility','off');
        plot(ax, ox(end), oy(end), 'k+','MarkerSize',6,'LineWidth',1.0,'HandleVisibility','off');
    end

    plot(ax, nan, nan, 'k--','LineWidth',1.0,'DisplayName','Obstacle center path');
    plot(ax, nan, nan, 'k:','LineWidth',1.2,'DisplayName','Keep-out at $t=0$');
    plot(ax, nan, nan, 'k-','LineWidth',1.2,'DisplayName','Keep-out at $t=T$');
end

function draw_circle(ax, cx, cy, R, lineStyle, lineWidth)
    rectangle(ax, 'Position',[cx-R, cy-R, 2*R, 2*R], ...
        'Curvature',[1 1], ...
        'EdgeColor','k', ...
        'LineStyle', lineStyle, ...
        'LineWidth', lineWidth);
end

function [minClear, violMask, xD, yD] = dense_clearance_check(solution, problem, nSub)
    p  = problem.params;
    tN = solution.time_vec(:)';
    xN = solution.states(1,:); yN = solution.states(2,:);
    vxN = solution.states(3,:); vyN = solution.states(4,:);

    if ~isfield(p,'obs') || isempty(p.obs)
        minClear = inf; violMask = false(1,numel(tN));
        xD = xN; yD = yN;
        return;
    end

    N = numel(tN)-1;
    xD = []; yD = []; tD = [];

    for k = 1:N
        s = linspace(0,1,nSub+1);
        if k < N, s = s(1:end-1); end

        tk  = tN(k);
        tk1 = tN(k+1);
        h = (tk1 - tk);

        t_local = tk + (tk1 - tk)*s;

        [x_local, ~] = hermite_1d(xN(k),  vxN(k),  xN(k+1),  vxN(k+1),  h, s);
        [y_local, ~] = hermite_1d(yN(k),  vyN(k),  yN(k+1),  vyN(k+1),  h, s);

        tD = [tD, t_local]; %#ok<AGROW>
        xD = [xD, x_local]; %#ok<AGROW>
        yD = [yD, y_local]; %#ok<AGROW>
    end

    clear_all = inf(size(tD));
    for i = 1:numel(p.obs)
        ox = p.obs(i).pos0(1) + p.obs(i).vel(1)*tD;
        oy = p.obs(i).pos0(2) + p.obs(i).vel(2)*tD;
        dist = sqrt((xD-ox).^2 + (yD-oy).^2);
        clear_i = dist - (p.obs(i).R + p.safety_margin);
        clear_all = min(clear_all, clear_i);
    end

    minClear = min(clear_all);
    violMask = (clear_all < 0);
end

function [p_s, v_s] = hermite_1d(p0, v0, p1, v1, h, s)
    s = s(:).'; % 1xM
    h00 =  2*s.^3 - 3*s.^2 + 1;
    h10 =      s.^3 - 2*s.^2 + s;
    h01 = -2*s.^3 + 3*s.^2;
    h11 =      s.^3 - s.^2;

    p_s = p0*h00 + (h*v0)*h10 + p1*h01 + (h*v1)*h11;

    dh00 = 6*s.^2 - 6*s;
    dh10 = 3*s.^2 - 4*s + 1;
    dh01 = -dh00;
    dh11 = 3*s.^2 - 2*s;

    dp_ds = p0*dh00 + (h*v0)*dh10 + p1*dh01 + (h*v1)*dh11;
    v_s = (1/h) * dp_ds;
end

function [xS,yS,xG,yG] = get_start_goal(solution, problem)
    x = solution.states(1,:); y = solution.states(2,:);
    if isfield(problem,'bounds') && isfield(problem.bounds,'initial_states')
        xS = problem.bounds.initial_states(1); yS = problem.bounds.initial_states(2);
    else
        xS = x(1); yS = y(1);
    end
    if isfield(problem,'bounds') && isfield(problem.bounds,'terminal_states') && ...
            ~any(isnan(problem.bounds.terminal_states(1:2)))
        xG = problem.bounds.terminal_states(1); yG = problem.bounds.terminal_states(2);
    else
        xG = x(end); yG = y(end);
    end
end

function [maxIneq, maxEq, constrV] = extract_solver_feas_text(summary)
    % New solver fields
    if isfield(summary,'max_node_ineq_violation')
        maxIneq = summary.max_node_ineq_violation;
    elseif isfield(summary,'max_ineq_violation')
        maxIneq = summary.max_ineq_violation;
    else
        maxIneq = NaN;
    end

    if isfield(summary,'max_node_eq_violation')
        maxEq = summary.max_node_eq_violation;
    elseif isfield(summary,'max_eq_violation')
        maxEq = summary.max_eq_violation;
    else
        maxEq = NaN;
    end

    if isfield(summary,'constrviolation')
        constrV = summary.constrviolation;
    else
        constrV = max([maxIneq, maxEq]);
    end
end

function s = set_defaults(s, d)
    f = fieldnames(d);
    for k = 1:numel(f)
        if ~isfield(s,f{k}) || isempty(s.(f{k}))
            s.(f{k}) = d.(f{k});
        end
    end
end

function v = getfield_default(s, f, d)
    if isstruct(s) && isfield(s,f) && ~isempty(s.(f)), v = s.(f);
    else, v = d; end
end
