clc; clear; close all
%% Parameters
Ls = [0.1, 0.15, 0.1, 0.15, 0.1]; % section lengths [m]
m = numel(Ls);
% Target tip position [X Y Z]
r_des = [0.3 0.15 0.15];
% Obstacle points to stay away from
r_obs1 = [0.242 -0.0459 -0.021];
d_safe1 = 0.05;
r_obs2 = [0.15 0.20 0.10];     % ← example second obstacle, change as needed
d_safe2 = 0.025;               % ← example second distance, change as needed
% Material & cross-section properties
E_mod = 210e9; % Young's modulus [Pa]
I_area = 1e-8; % area moment of inertia [m⁴]
nu = 0.3;
G_mod = E_mod / (2*(1 + nu));
J_area = 2 * I_area; % torsional constant [m⁴]
%% Initial guess & bounds
kappa0 = 1.0 * ones(1, m); % initial curvature guess [1/m] — lowered a bit
phi0 = zeros(1, m); % initial twist angles [rad]
lb_k = -12; ub_k = 12; % wider curvature bounds — helps reachability
lb_p = -pi; ub_p = pi;
%% Run optimization
[kOpt, phiOpt, Emin, tip_err] = optimiseMinimalEnergyRod( ...
    Ls, r_des, r_obs1, d_safe1, r_obs2, d_safe2, kappa0, phi0, ...
    E_mod, I_area, G_mod, J_area, ...
    lb_k, ub_k, lb_p, ub_p, []);
fprintf('Tip position error = %.2e m\n', tip_err);
fprintf('Elastic energy = %.4f J\n', Emin);
%% Plot final shape
figure('Color','w'); hold on;
% Target tip
plot3(r_des(2), r_des(3), r_des(1), '*r', 'MarkerSize',12, 'LineWidth',2.5, 'DisplayName','Target tip');
hold on;
plot3(r_obs1(2), r_obs1(3), r_obs1(1), '*k', 'MarkerSize', 8, 'LineWidth', 1.5);
plot3(r_obs2(2), r_obs2(3), r_obs2(1), 'xk', 'MarkerSize', 8, 'LineWidth', 1.5);
% Optimized rod
nPlot = ceil(400 * Ls / max(Ls)); % denser sampling for smooth plot
[~, rLine, r_junc] = rodFK(kOpt, phiOpt, Ls, nPlot);
plot3(r_junc(:,2), r_junc(:,3), r_junc(:,1), '.g', 'MarkerSize',24, 'LineWidth',4, 'DisplayName','Junctions');
plot3(rLine(:,2), rLine(:,3), rLine(:,1), 'r-', 'LineWidth',3, 'DisplayName','Centerline');
grid on; axis equal; box on;
xlabel('Y (m)'); ylabel('Z (m)'); zlabel('X (m)');
set(gca, 'FontSize',13);
xlim([-0.6 0.6]); ylim([-0.6 0.6]); zlim([0 0.8]);
legend('Location','bestoutside');
title('Minimum Energy Shape Reaching Target Tip','Interpreter','latex','FontSize',15);
hold off;
%% =======================================================================
% MAIN OPTIMIZER – only tip position constraint
% =======================================================================
function [kappa_opt, phi_opt, E_opt, tip_err] = optimiseMinimalEnergyRod( ...
    Ls, r_des,r_obs1,d_safe1,r_obs2,d_safe2, kappa0, phi0, ...
    E_mod, I_area, G_mod, J_area, ...
    lb_kappa, ub_kappa, lb_phi, ub_phi,nPts_vec)
% ---------- sizes & default sampling ----------
Ls = Ls(:)'; m = numel(Ls);
if nargin<13 || isempty(nPts_vec)
    nPts_vec = max(50, ceil(200*Ls/max(Ls)));
end
% ---------- variable pack ----------
x0 = [kappa0(:).' , phi0(:).']; % 1×(2m-1)
lb = [lb_kappa(:).', lb_phi(:).'];
ub = [ub_kappa(:).', ub_phi(:).'];
% Objective = elastic energy only
    obj = @(x) energyObjective(x, Ls, E_mod, I_area, G_mod, J_area);
% Constraint: tip position equality + two clearances
    nonlcon = @(x) tipConstraint(x, Ls, nPts_vec, r_des,r_obs1,d_safe1,r_obs2,d_safe2);
opts = optimoptions('fmincon','Algorithm','sqp', ...
'Display','iter','MaxFunctionEvaluations',1e4);
    [x_opt, E_opt] = fmincon(obj, x0, [], [], [], [], lb, ub, nonlcon, opts);
    kappa_opt = x_opt(1:m);
    phi_opt = x_opt(m+1:end);
% Compute tip error for reporting
    [~, r_line] = rodFK(kappa_opt, phi_opt, Ls, nPts_vec);
    tip_err = norm(r_line(end,:) - r_des(:)');
end
%% =======================================================================
% ELASTIC ENERGY (objective)
% =======================================================================
function E = energyObjective(x, Ls, E_mod, I_area, G_mod, J_area)
    m = numel(Ls);
    kappa = x(1:m);
    phi = x(m+1:end);
    E_bend = 0.5 * E_mod * I_area * sum(Ls .* kappa.^2);
    E_twist = 0.5 * G_mod * J_area * sum(phi.^2);
    E = E_bend + E_twist;
end
%% =======================================================================
% CONSTRAINT: tip equality + two clearances
% =======================================================================
function [c, ceq] = tipConstraint(x, Ls, nPts_vec, r_des,r_obs1,d_safe1,r_obs2,d_safe2)
    m = numel(Ls);
    kappa = x(1:m);
    phi = x(m+1:end);
    [~, r_line] = rodFK(kappa, phi, Ls, nPts_vec);
    ceq = r_line(end,:).' - r_des(:); % equality constraint
    % Inequalities for two obstacles
    distances1 = vecnorm(r_line - r_obs1(:)', 2, 2);
    dmin1 = min(distances1);
    distances2 = vecnorm(r_line - r_obs2(:)', 2, 2);
    dmin2 = min(distances2);
    c = [d_safe1 - dmin1; 
         d_safe2 - dmin2]; % both <= 0 means both clearances satisfied
end
%% =======================================================================
% FORWARD KINEMATICS (unchanged)
% =======================================================================
function [s_tot, r_tot, r_junc] = rodFK(kappa, phi, Ls, nPts_vec)
    m = numel(kappa);
if numel(phi) ~= m || numel(Ls) ~= m || numel(nPts_vec) ~= m
        error('Dimension mismatch in rodFK');
end
    curvature_funs = cellfun(@(k) @(s) k, num2cell(kappa), 'UniformOutput', false);
    torsion_funs = repmat({@(s) 0}, 1, m);
    T0 = [1 0 0];
    Nor0 = [0 1 0]*cos(phi(1)) + [0 0 1]*sin(phi(1));
    B0 = cross(T0, Nor0);
    [s_tot, r_tot, ~, ~, ~] = rod_multi_sections_customStart( ...
        curvature_funs, torsion_funs, Ls, nPts_vec, phi(2:end), ...
        [0 0 0], T0, Nor0, B0);
    cumL = [0 cumsum(Ls)];
    r_junc = zeros(numel(cumL), 3);
for j = 1:numel(cumL)
        [~, idx] = min(abs(s_tot - cumL(j)));
        r_junc(j,:) = r_tot(idx,:);
end
end
function [s_tot, r_tot, T_tot, Nor_tot, B_tot] = rod_multi_sections_customStart( ...
    curvature_funs, torsion_funs, Ls, nPts_vec, phis, r0, T0, N0, B0)
    m = numel(Ls);
    s_tot = []; r_tot = []; T_tot = []; Nor_tot = []; B_tot = [];
    r_prev = r0; T_prev = T0; Nor_prev = N0; B_prev = B0; cumL = 0;
for i = 1:m
        Li = Ls(i); ni = nPts_vec(i);
        [s_i, r_i, T_i, Nor_i, B_i] = rod_segment( ...
            curvature_funs{i}, torsion_funs{i}, Li, ni, r_prev, T_prev, Nor_prev, B_prev);
        s_i_shift = cumL + s_i;
if i == 1
            s_tot = s_i_shift; r_tot = r_i; T_tot = T_i; Nor_tot = Nor_i; B_tot = B_i;
else
            s_tot = [s_tot; s_i_shift(2:end)];
            r_tot = [r_tot; r_i(2:end,:)];
            T_tot = [T_tot; T_i(2:end,:)];
            Nor_tot = [Nor_tot; Nor_i(2:end,:)];
            B_tot = [B_tot; B_i(2:end,:)];
end
        cumL = cumL + Li;
if i < m
            phi = phis(i);
            T_prev = T_i(end,:);
            Nend = Nor_i(end,:);
            Bend = B_i(end,:);
            Nor_prev = cos(phi)*Nend + sin(phi)*Bend;
            B_prev = -sin(phi)*Nend + cos(phi)*Bend;
            r_prev = r_i(end,:);
end
end
end
function [s_vec, r, T, Nor, B] = rod_segment(curv_fun, tau_fun, L_seg, nPts, r0, T0, N0, B0)
    s_vec = linspace(0, L_seg, nPts).';
    r = zeros(nPts,3); r(1,:) = r0;
    T = zeros(nPts,3); T(1,:) = T0;
    Nor = zeros(nPts,3); Nor(1,:) = N0;
    B = zeros(nPts,3); B(1,:) = B0;
for k = 1:nPts-1
        ds = s_vec(k+1) - s_vec(k);
        s = s_vec(k);
        kappa = curv_fun(s);
        tau = tau_fun(s);
        y0 = [T(k,:)'; Nor(k,:)'; B(k,:)'];
        f = @(y) [kappa*y(4:6); -kappa*y(1:3) + tau*y(7:9); -tau*y(4:6)];
        k1 = f(y0);
        k2 = f(y0 + 0.5*ds*k1);
        k3 = f(y0 + 0.5*ds*k2);
        k4 = f(y0 + ds*k3);
        y1 = y0 + (ds/6)*(k1 + 2*k2 + 2*k3 + k4);
        T(k+1,:) = y1(1:3)' / norm(y1(1:3));
        Nor(k+1,:) = y1(4:6)' / norm(y1(4:6));
        B(k+1,:) = cross(T(k+1,:), Nor(k+1,:));
        B(k+1,:) = B(k+1,:) / norm(B(k+1,:));
        Nor(k+1,:) = cross(B(k+1,:), T(k+1,:));
        r(k+1,:) = r(k,:) + ds * T(k,:);
end
end