% =========================================================================
%  KIRCHHOFF ROD THROUGH A CYLINDER — TIP POSITION TARGETING
%  -----------------------------------------------------------------
%  Theory  : Kirchhoff / Cosserat rod with Darboux frame kinematics
%  Solver  : Augmented Lagrangian (tip) + log-barrier (wall) + fminunc
%  Outputs : curvature κ(s), torsion τ(s), 3-D centreline r(s)
%
%  Variables per segment i  :  [k1_i, k2_i, tau_i]
%     k1, k2  — material bending curvatures  (1/m)
%     tau     — material twist / torsion      (1/m)
%
%  State vector at each s   :  [r(3); d1(3); d2(3); d3(3)]  (12 DOF)
%
%  Usage:
%     kirchhoff_rod_cylinder          % runs with built-in defaults
%     kirchhoff_rod_cylinder(params)  % pass a struct to override anything
%
%  Tunable parameters are collected in the PARAMETERS section below.
% =========================================================================
clc; close all;
clear;

%% ── PARAMETERS ──────────────────────────────────────────────────────────
p.N          = 60;               % number of rod segments
p.L          = 0.5;              % total rod length  [m]

% Material (isotropic circular cross-section)
p.EI         = 1e-3;             % bending stiffness  [N·m²]
p.GJ         = 5e-4;             % torsional stiffness [N·m²]

% Cylinder
p.R_cyl      = 0.08;            % cylinder inner radius  [m]
p.L_cyl      = 0.2;             % cylinder length  [m]  (must be < L)
p.axis_start = [0.0; 0; 0.0];     % point on cylinder axis at entry  [m]
p.axis_dir   = [0; 0; 1];        % unit vector along cylinder axis

% Desired tip position
p.r_tip_des  = [0.05; 0.02; 0.35]; % [m]

% Initial base conditions
p.r0         = [0; 0; 0];        % base position
p.d10        = [1; 0; 0];        % initial d1 (material frame)
p.d20        = [0; 1; 0];        % initial d2
p.d30        = [0; 0; 1];        % initial tangent d3 = r'(0)

% Barrier / penalty weights
p.mu_wall    = 1e6;              % log-barrier weight for wall avoidance
p.mu_tip     = 1e5;              % tip penalty (augmented Lagrangian)
p.lambda0    = zeros(3,1);       % initial Lagrange multipliers (tip)
p.n_AL       = 8;                % augmented Lagrangian outer iterations
p.mu_scale   = 5;                % factor to increase mu_tip each AL step

% Optimiser
p.max_iter   = 400;
p.tol_fun    = 1e-10;
p.tol_x      = 1e-10;
p.verbose    = true;

kirchhoff_rod_cylinder(p)       % pass your struct in
function kirchhoff_rod_cylinder(user_params)


% Override with user-supplied struct if provided
if nargin > 0
    fnames = fieldnames(user_params);
    for k = 1:numel(fnames)
        p.(fnames{k}) = user_params.(fnames{k});
    end
end

% Derived quantities
p.ds       = p.L / p.N;
p.s        = linspace(0, p.L, p.N+1);
p.axis_dir = p.axis_dir / norm(p.axis_dir);

%% ── INITIAL GUESS ────────────────────────────────────────────────────────
x0 = zeros(3 * p.N, 1);   % [k1_1..k1_N, k2_1..k2_N, tau_1..tau_N]

%% ── AUGMENTED LAGRANGIAN OUTER LOOP ──────────────────────────────────────
lambda = p.lambda0;
mu_tip = p.mu_tip;
x      = x0;

opts = optimoptions('fminunc', ...
    'Algorithm',                'quasi-newton', ...
    'SpecifyObjectiveGradient', false, ...
    'MaxIterations',            p.max_iter, ...
    'FunctionTolerance',        p.tol_fun, ...
    'StepTolerance',            p.tol_x, ...
    'Display',                  'off');

if p.verbose
    fprintf('\n%-5s  %-14s  %-14s  %-14s\n', ...
        'AL it', 'E_elastic', 'tip_error', '|lambda|');
    fprintf('%s\n', repmat('-',1,55));
end

for al_iter = 1:p.n_AL

    cost_fn  = @(xv) total_cost(xv, lambda, mu_tip, p);
    [x, ~]   = fminunc(cost_fn, x, opts);

    % Evaluate tip residual after this AL step
    [r_all, ~] = shoot_rod(x, p);
    delta_r    = r_all(:,end) - p.r_tip_des;

    % Multiplier + penalty update
    lambda = lambda + mu_tip * delta_r;
    mu_tip = mu_tip * p.mu_scale;

    E_el = elastic_energy(x, p);
    if p.verbose
        fprintf('%-5d  %-14.6e  %-14.6e  %-14.6e\n', ...
            al_iter, E_el, norm(delta_r), norm(lambda));
    end

    if norm(delta_r) < 1e-6
        if p.verbose
            fprintf('Converged: tip error = %.2e m\n', norm(delta_r));
        end
        break
    end
end

%% ── POST-PROCESS & VISUALISE ─────────────────────────────────────────────
[r_all, frames] = shoot_rod(x, p);

k1    = x(1:p.N);
k2    = x(p.N+1:2*p.N);
tau   = x(2*p.N+1:3*p.N);
kappa = sqrt(k1.^2 + k2.^2);

s_mid = (p.s(1:end-1) + p.s(2:end)) / 2;

% Figure 1 — 3D shape
figure('Name','Kirchhoff Rod - 3D Shape','NumberTitle','off','Color','w');
ax = axes; hold on; grid on; axis equal;
xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');
title('Kirchhoff rod through cylinder to tip target','FontWeight','normal');

draw_cylinder(p, ax);
plot3(r_all(1,:), r_all(2,:), r_all(3,:), 'b-', 'LineWidth', 2.5, ...
    'DisplayName','rod centreline');
scatter3(p.r0(1), p.r0(2), p.r0(3), 80, 'k', 'filled', ...
    'DisplayName','base');
scatter3(p.r_tip_des(1), p.r_tip_des(2), p.r_tip_des(3), 120, 'r', ...
    'p', 'filled', 'DisplayName','desired tip r*');
scatter3(r_all(1,end), r_all(2,end), r_all(3,end), 80, 'g', ...
    'filled', 'DisplayName','achieved tip');

skip = max(1, floor(p.N/10));
for i = 1:skip:p.N
    r_i = r_all(:, i);
    sc  = p.ds * 3;
    quiver3(r_i(1),r_i(2),r_i(3), ...
        frames(1,i)*sc, frames(2,i)*sc, frames(3,i)*sc, ...
        0, 'r', 'LineWidth', 1.2, 'HandleVisibility', 'off');
    quiver3(r_i(1),r_i(2),r_i(3), ...
        frames(4,i)*sc, frames(5,i)*sc, frames(6,i)*sc, ...
        0, 'g', 'LineWidth', 1.2, 'HandleVisibility', 'off');
    quiver3(r_i(1),r_i(2),r_i(3), ...
        frames(7,i)*sc, frames(8,i)*sc, frames(9,i)*sc, ...
        0, 'b', 'LineWidth', 1.2, 'HandleVisibility', 'off');
end
legend('Location','best'); view(45,25);

% Figure 2 — strain profiles
figure('Name','Kirchhoff Rod - Strain Profiles','NumberTitle','off','Color','w');
subplot(3,1,1);
plot(s_mid, k1, 'b-o', 'MarkerSize', 4); grid on;
xlabel('arc length s [m]'); ylabel('\kappa_1 [1/m]');
title('Bending curvature \kappa_1(s)','FontWeight','normal');

subplot(3,1,2);
plot(s_mid, k2, 'r-o', 'MarkerSize', 4); grid on;
xlabel('arc length s [m]'); ylabel('\kappa_2 [1/m]');
title('Bending curvature \kappa_2(s)','FontWeight','normal');

subplot(3,1,3);
plot(s_mid, tau, 'm-o', 'MarkerSize', 4); grid on;
xlabel('arc length s [m]'); ylabel('\tau [1/m]');
title('Torsion / twist \tau(s)','FontWeight','normal');

% Figure 3 — total curvature
figure('Name','Kirchhoff Rod - Total Curvature','NumberTitle','off','Color','w');
plot(s_mid, kappa, 'k-', 'LineWidth', 2); grid on;
xlabel('arc length s [m]');
ylabel('\kappa(s) [1/m]');
title('Total curvature along rod','FontWeight','normal');

% Summary
fprintf('\n--- Solution summary -------------------------------------------\n');
fprintf('  Elastic energy         : %.6e N.m\n', ...
    sum(0.5*p.ds*(p.EI*(k1.^2+k2.^2) + p.GJ*tau.^2)));
fprintf('  Max total curvature    : %.4f 1/m\n', max(kappa));
fprintf('  Max torsion |tau|      : %.4f 1/m\n', max(abs(tau)));
fprintf('  Tip achieved  [m]      : [%.4f, %.4f, %.4f]\n', ...
    r_all(1,end), r_all(2,end), r_all(3,end));
fprintf('  Tip desired   [m]      : [%.4f, %.4f, %.4f]\n', p.r_tip_des');
fprintf('  Tip error              : %.2e m\n', norm(r_all(:,end)-p.r_tip_des));
fprintf('----------------------------------------------------------------\n\n');

end  % ── end main function ──────────────────────────────────────────────


%% =========================================================================
%  TOTAL COST  (objective for fminunc)
% =========================================================================
function [F, G] = total_cost(x, lambda, mu_tip, p)
    E_el       = elastic_energy(x, p);
    E_wall     = wall_barrier(x, p);
    [r_all, ~] = shoot_rod(x, p);
    delta_r    = r_all(:,end) - p.r_tip_des;

    E_tip = lambda' * delta_r + 0.5 * mu_tip * (delta_r' * delta_r);
    F     = E_el + E_wall + E_tip;

    if nargout > 1
        G = numerical_gradient(x, lambda, mu_tip, p);
    end
end


%% =========================================================================
%  ELASTIC ENERGY
% =========================================================================
function E = elastic_energy(x, p)
    N   = p.N;
    k1  = x(1:N);
    k2  = x(N+1:2*N);
    tau = x(2*N+1:3*N);
    E   = 0.5 * p.ds * sum(p.EI*(k1.^2 + k2.^2) + p.GJ*tau.^2);
end


%% =========================================================================
%  LOG-BARRIER WALL AVOIDANCE
% =========================================================================
function E_wall = wall_barrier(x, p)
    [r_all, ~] = shoot_rod(x, p);
    E_wall = 0;

    for i = 1:p.N+1
        proj  = (r_all(:,i) - p.axis_start)' * p.axis_dir;
        if proj >= 0 && proj <= p.L_cyl
            r_vec = r_all(:,i) - p.axis_start - proj * p.axis_dir;
            rho   = norm(r_vec);
            gap   = p.R_cyl - rho;
            if gap <= 0
                E_wall = E_wall + p.mu_wall * gap^2;
            else
                E_wall = E_wall - p.mu_wall * p.ds * log(gap / p.R_cyl);
            end
        end
    end
end


%% =========================================================================
%  FORWARD SHOOT — RK4 integration of Darboux ODE
%  r_all  : 3 x (N+1)   centreline positions
%  frames : 9 x N       material frame [d1;d2;d3] at each segment start
% =========================================================================
function [r_all, frames] = shoot_rod(x, p)
    N  = p.N;
    ds = p.ds;

    k1  = x(1:N);
    k2  = x(N+1:2*N);
    tau = x(2*N+1:3*N);

    state = [p.r0; p.d10; p.d20; p.d30];   % 12 x 1

    r_all  = zeros(3, N+1);
    frames = zeros(9, N);
    r_all(:,1) = state(1:3);

    for i = 1:N
        kv  = [k1(i); k2(i); tau(i)];

        k_1 = darboux_rhs(state,               kv);
        k_2 = darboux_rhs(state + 0.5*ds*k_1,  kv);
        k_3 = darboux_rhs(state + 0.5*ds*k_2,  kv);
        k_4 = darboux_rhs(state +     ds*k_3,  kv);

        state = state + (ds/6)*(k_1 + 2*k_2 + 2*k_3 + k_4);
        state = orthonormalise_frame(state);

        r_all(:,i+1)  = state(1:3);
        frames(1:3,i) = state(4:6);    % d1
        frames(4:6,i) = state(7:9);    % d2
        frames(7:9,i) = state(10:12);  % d3
    end
end


%% =========================================================================
%  DARBOUX ODE RHS
%  r'  = d3
%  d1' = tau*d2 - k2*d3
%  d2' = k1*d3  - tau*d1
%  d3' = k2*d1  - k1*d2
% =========================================================================
function drhs = darboux_rhs(state, kv)
    d1  = state(4:6);
    d2  = state(7:9);
    d3  = state(10:12);
    k1  = kv(1);  k2 = kv(2);  tau = kv(3);

    drhs = [d3; ...
            tau*d2 - k2*d3; ...
            k1*d3  - tau*d1; ...
            k2*d1  - k1*d2];
end


%% =========================================================================
%  GRAM-SCHMIDT ORTHONORMALISATION (prevents frame drift)
% =========================================================================
function state = orthonormalise_frame(state)
    r  = state(1:3);
    d1 = state(4:6);
    d2 = state(7:9);
    d3 = state(10:12);

    d3 = d3 / norm(d3);
    d1 = d1 - (d1'*d3)*d3;
    d1 = d1 / norm(d1);
    d2 = cross(d3, d1);
    d2 = d2 / norm(d2);

    state = [r; d1; d2; d3];
end


%% =========================================================================
%  NUMERICAL GRADIENT (central finite differences)
% =========================================================================
function G = numerical_gradient(x, lambda, mu_tip, p)
    n = numel(x);
    G = zeros(n, 1);
    h = 1e-6;

    for i = 1:n
        xp = x;  xp(i) = xp(i) + h;
        xm = x;  xm(i) = xm(i) - h;

        [rp, ~] = shoot_rod(xp, p);
        [rm, ~] = shoot_rod(xm, p);

        drp = rp(:,end) - p.r_tip_des;
        drm = rm(:,end) - p.r_tip_des;

        Fp = elastic_energy(xp,p) + wall_barrier(xp,p) ...
           + lambda'*drp + 0.5*mu_tip*(drp'*drp);

        Fm = elastic_energy(xm,p) + wall_barrier(xm,p) ...
           + lambda'*drm + 0.5*mu_tip*(drm'*drm);

        G(i) = (Fp - Fm) / (2*h);
    end
end


%% =========================================================================
%  DRAW CYLINDER (transparent surface mesh)
% =========================================================================
function draw_cylinder(p, ax)
    theta  = linspace(0, 2*pi, 40);
    t_axis = linspace(0, p.L_cyl, 20);

    ax_dir = p.axis_dir;
    if abs(ax_dir(1)) < 0.9
        perp1 = cross(ax_dir, [1;0;0]);
    else
        perp1 = cross(ax_dir, [0;1;0]);
    end
    perp1 = perp1 / norm(perp1);
    perp2 = cross(ax_dir, perp1);

    X = zeros(numel(theta), numel(t_axis));
    Y = X;  Z = X;

    for ti = 1:numel(t_axis)
        ctr = p.axis_start + t_axis(ti)*ax_dir;
        for tj = 1:numel(theta)
            pt = ctr + p.R_cyl*(cos(theta(tj))*perp1 + sin(theta(tj))*perp2);
            X(tj,ti) = pt(1);
            Y(tj,ti) = pt(2);
            Z(tj,ti) = pt(3);
        end
    end

    surf(ax, X, Y, Z, ...
        'FaceAlpha', 0.15, 'FaceColor', [0.8 0.85 0.9], ...
        'EdgeColor', [0.6 0.6 0.6], 'EdgeAlpha', 0.3, ...
        'DisplayName', 'cylinder wall');

    for cap = [0, p.L_cyl]
        ctr = p.axis_start + cap*ax_dir;
        xc  = zeros(1, numel(theta));
        yc  = xc;  zc = xc;
        for tj = 1:numel(theta)
            pt = ctr + p.R_cyl*(cos(theta(tj))*perp1 + sin(theta(tj))*perp2);
            xc(tj) = pt(1);  yc(tj) = pt(2);  zc(tj) = pt(3);
        end
        plot3(ax, xc, yc, zc, 'Color', [0.5 0.5 0.5], ...
            'LineWidth', 1, 'HandleVisibility', 'off');
    end
end
