clc;
clear;
close all;

%% ========================================================================
%  TDCR PCC OPTIMISATION  —  SPHERICAL OBSTACLE AVOIDANCE
%  ---------------------------------------------------------------
%  METHOD:
%
%   TDCR Piecewise Constant Curvature
% + Frenet Frame Integration
% + Geometric RK4 Integration
% + Log Barrier Sphere Avoidance   <-- KEY CHANGE
%     Feasible region : OUTSIDE the sphere  (rho >= R + d_safe)
%     Barrier blows up as the rod approaches the sphere surface
% + Augmented Lagrangian Tip Targeting
%
%  Unknowns:
%       x = [kappa_1 ... kappa_m  phi_1 ... phi_m]
%
%  Optimisation:
%
%     F(x) =
%       E_elastic
%     + E_sphere           <-- log barrier: penalises entry into sphere
%     + E_tip(AL)
%
%  Uses:
%       fminunc (unconstrained optimisation)
%
%  GEOMETRY NOTE:
%     The rod must navigate AROUND the sphere to reach the target.
%     The sphere sits between the robot base and the target, so the
%     optimiser is forced to curve the backbone around the obstacle.
%
%% ========================================================================

%% ========================================================================
%  TDCR PARAMETERS
%% ========================================================================

Ls = [0.10 0.10 0.10 0.10 0.10];   % section arc-lengths [m]

m = numel(Ls);

Ltot = sum(Ls);

%% ========================================================================
%  TARGET TIP
%% ========================================================================
%  Place target beyond / beside the sphere so the robot must go around it.

r_des = [0.1  0.01  0.30];

%% ========================================================================
%  MATERIAL
%% ========================================================================

EI = 1.0;    % bending stiffness  [N·m²]
GJ = 0.1;    % torsional stiffness [N·m²]

%% ========================================================================
%  SPHERICAL OBSTACLE
%% ========================================================================
%
%  Feasible set:  { p ∈ R³  |  ‖p − c‖ ≥ R + d_safe }
%  i.e. the robot must stay OUTSIDE the sphere.
%
%  The log-barrier is:
%
%     Φ(gap)  =  −μ · ln( gap / R )        if  gap > 0   (feasible)
%             =   10⁶ · gap²               if  gap ≤ 0   (collision)
%
%  where  gap  =  ‖p − c‖ − R − d_safe   (positive when outside)
%

sph.center  = [0.02  0.01  0.20];   % sphere centre [m]
sph.radius  = 0.06;                 % sphere radius [m]

d_safe = 0.004;                     % clearance buffer [m]

%% ========================================================================
%  BARRIER + AL PARAMETERS
%% ========================================================================

mu_sphere  = 5e-3;    % log-barrier weight for sphere

mu_tip     = 10;      % initial AL penalty weight

mu_scale   = 5;       % AL penalty growth factor per outer iteration

lambda     = zeros(3,1);   % Lagrange multiplier (3-D tip error)

nAL        = 10;           % max augmented Lagrangian outer iterations

%% ========================================================================
%  DISCRETISATION
%% ========================================================================

nPts = 80 * ones(1,m);   % integration points per section

%% ========================================================================
%  INITIAL GUESS
%% ========================================================================
%  Small non-zero curvatures with a slight bending-plane offset help the
%  optimiser find a path that curves away from the sphere immediately.

x0 = [1.5*ones(1,m)   0.3*ones(1,m)];

x  = x0;

%% ========================================================================
%  OPTIMISER OPTIONS
%% ========================================================================

opts = optimoptions('fminunc', ...
    'Algorithm',              'quasi-newton', ...
    'Display',                'iter',         ...
    'MaxIterations',          600,            ...
    'MaxFunctionEvaluations', 3e5,            ...
    'OptimalityTolerance',    1e-8,           ...
    'StepTolerance',          1e-10);

%% ========================================================================
%  AUGMENTED LAGRANGIAN OUTER LOOP
%% ========================================================================

fprintf('\n');
fprintf('=====================================================\n');
fprintf(' TDCR PCC + Sphere Barrier + Augmented Lagrangian\n');
fprintf('=====================================================\n\n');

for al = 1:nAL

    fprintf('AL ITERATION %d\n', al);

    obj = @(xx) totalCost( ...
        xx,        ...
        Ls,        ...
        nPts,      ...
        EI,        ...
        GJ,        ...
        sph,       ...
        d_safe,    ...
        r_des,     ...
        mu_sphere, ...
        lambda,    ...
        mu_tip);

    x = fminunc(obj, x, opts);

    %% ---- evaluate tip error at current solution

    [~, rline] = rodFK(x(1:m), x(m+1:end), Ls, nPts);

    e = rline(end,:)' - r_des(:);

    fprintf('Tip Error = %.6e\n', norm(e));

    %% ---- Augmented Lagrangian dual + penalty update

    lambda = lambda + mu_tip * e;

    mu_tip = mu_tip * mu_scale;

    fprintf('mu_tip    = %.3e\n', mu_tip);
    fprintf('--------------------------------------------\n');

    if norm(e) < 1e-5
        fprintf('Converged.\n');
        break;
    end
end

%% ========================================================================
%  FINAL RESULT
%% ========================================================================

kappa = x(1:m);
phi   = x(m+1:end);

[s, rline, rj] = rodFK(kappa, phi, Ls, nPts);

Eel = elasticEnergy(kappa, phi, Ls, EI, GJ);

%% ---- minimum clearance from sphere surface

dists = vecnorm(rline - sph.center, 2, 2);   % distance of each point to centre
min_clearance = min(dists) - sph.radius;

fprintf('\n');
fprintf('=====================================================\n');
fprintf('FINAL RESULT\n');
fprintf('=====================================================\n');
fprintf('Elastic Energy    = %.6f\n', Eel);
fprintf('Tip Error         = %.6e\n', norm(rline(end,:) - r_des));
fprintf('Min Sphere Gap    = %.4f m  (positive = no collision)\n', min_clearance);
fprintf('\n');

%% ========================================================================
%  VISUALISATION
%% ========================================================================

figure('Color', 'w');
hold on;
grid on;
axis equal;
view(3);

xlabel('X [m]');
ylabel('Y [m]');
zlabel('Z [m]');

title('TDCR PCC Optimisation — Sphere Avoidance (Outside Feasible)');

%% ---- sphere surface

drawSphere(sph.center, sph.radius);

%% ---- safety-margin sphere (wireframe)

drawSphereWire(sph.center, sph.radius + d_safe);

%% ---- robot backbone

plot3(rline(:,1), rline(:,2), rline(:,3), ...
    'r-', 'LineWidth', 3);

%% ---- section junctions

plot3(rj(:,1), rj(:,2), rj(:,3), ...
    '.k', 'MarkerSize', 22);

%% ---- target point

plot3(r_des(1), r_des(2), r_des(3), ...
    'bp', 'MarkerSize', 16, 'MarkerFaceColor', 'b');

%% ---- sphere centre marker

plot3(sph.center(1), sph.center(2), sph.center(3), ...
    'k+', 'MarkerSize', 12, 'LineWidth', 2);

legend('Sphere (obstacle)', 'Safety margin', ...
       'Rod', 'Junctions', 'Target', 'Sphere centre', ...
       'Location', 'best');

%% ========================================================================
%  TOTAL COST  (objective passed to fminunc)
%% ========================================================================

function F = totalCost( ...
    x,        ...
    Ls,       ...
    nPts,     ...
    EI,       ...
    GJ,       ...
    sph,      ...
    d_safe,   ...
    r_des,    ...
    mu_sphere,...
    lambda,   ...
    mu_tip)

m     = numel(Ls);
kappa = x(1:m);
phi   = x(m+1:end);

%% ---- elastic energy

Eel = elasticEnergy(kappa, phi, Ls, EI, GJ);

%% ---- forward kinematics

[~, rline] = rodFK(kappa, phi, Ls, nPts);

%% ---- sphere barrier

Esph = sphereBarrier(rline, sph, d_safe, mu_sphere);

%% ---- tip augmented Lagrangian

e    = rline(end,:)' - r_des(:);
Etip = lambda' * e  +  0.5 * mu_tip * (e' * e);

%% ---- total cost

F = Eel + Esph + Etip;

fprintf('Eel=%7.4e  Esph=%7.4e  Etip=%7.4e\r', Eel, Esph, Etip);

end

%% ========================================================================
%  ELASTIC ENERGY
%% ========================================================================

function E = elasticEnergy(kappa, phi, Ls, EI, GJ)

Ebend  = 0.5 * EI * sum(Ls .* kappa.^2);
Etwist = 0.5 * GJ * sum(phi.^2);
E      = Ebend + Etwist;

end

%% ========================================================================
%  SPHERE LOG-BARRIER
%  ---------------------------------------------------------------
%  Feasible region:  OUTSIDE the sphere.
%
%    gap(p) = ‖p − c‖ − R − d_safe
%
%    gap > 0  →  outside (feasible)   →  −μ · ln( gap / R )
%    gap ≤ 0  →  inside  (collision)  →  10⁶ · gap²
%
%  The barrier → +∞ as the rod grazes the sphere surface (gap → 0⁺),
%  pushing the optimiser to keep the backbone outside.
%% ========================================================================

function Esph = sphereBarrier(rline, sph, d_safe, mu_sphere)

N  = size(rline, 1);
c  = sph.center(:)';
R  = sph.radius;

Esph = 0;

for i = 1:N

    p   = rline(i,:);
    rho = norm(p - c);          % distance from backbone point to sphere centre

    gap = rho - R - d_safe;     % positive when safely outside sphere

    if gap <= 0
        %% -- collision or too close: hard quadratic penalty
        Esph = Esph + 1e6 * gap^2;
    else
        %% -- log barrier: cost rises as gap → 0
        Esph = Esph - mu_sphere * log(gap / R);
    end

end

end

%% ========================================================================
%  FORWARD KINEMATICS
%% ========================================================================

function [s_tot, r_tot, r_junc] = rodFK(kappa, phi, Ls, nPts)

m = numel(kappa);

curv_funs = cellfun(@(k) @(s) k, num2cell(kappa), 'UniformOutput', false);
tors_funs = repmat({@(s) 0}, 1, m);

%% ---- initial frame aligned with z-axis, normal in bending plane of phi(1)

T0 = [0  0  1];
N0 = [1  0  0] * cos(phi(1)) + [0  1  0] * sin(phi(1));
B0 = cross(T0, N0);

[s_tot, r_tot, ~, ~, ~] = rodMulti( ...
    curv_funs, tors_funs, Ls, nPts, phi(2:end), ...
    [0 0 0], T0, N0, B0);

%% ---- extract junction positions

cumL   = [0  cumsum(Ls)];
r_junc = zeros(numel(cumL), 3);

for j = 1:numel(cumL)
    [~, idx]   = min(abs(s_tot - cumL(j)));
    r_junc(j,:) = r_tot(idx,:);
end

end

%% ========================================================================
%  MULTI-SECTION PROPAGATION
%% ========================================================================

function [s_tot, r_tot, T_tot, N_tot, B_tot] = rodMulti( ...
    curv_funs, tors_funs, Ls, nPts, phis, r0, T0, N0, B0)

m = numel(Ls);

s_tot = []; r_tot = []; T_tot = []; N_tot = []; B_tot = [];

r_prev = r0; T_prev = T0; N_prev = N0; B_prev = B0;
cumL   = 0;

for i = 1:m

    Li = Ls(i);
    ni = nPts(i);

    [s_i, r_i, T_i, N_i, B_i] = rodSegment( ...
        curv_funs{i}, tors_funs{i}, Li, ni, r_prev, T_prev, N_prev, B_prev);

    s_i = s_i + cumL;

    if i == 1
        s_tot = s_i; r_tot = r_i;
        T_tot = T_i; N_tot = N_i; B_tot = B_i;
    else
        s_tot = [s_tot;  s_i(2:end)];
        r_tot = [r_tot;  r_i(2:end,:)];
        T_tot = [T_tot;  T_i(2:end,:)];
        N_tot = [N_tot;  N_i(2:end,:)];
        B_tot = [B_tot;  B_i(2:end,:)];
    end

    cumL = cumL + Li;

    %% ---- rotate bending plane at junction

    if i < m
        ph     = phis(i);
        T_prev = T_i(end,:);
        Nend   = N_i(end,:);
        Bend   = B_i(end,:);

        N_prev =  cos(ph) * Nend + sin(ph) * Bend;
        B_prev = -sin(ph) * Nend + cos(ph) * Bend;
        r_prev =  r_i(end,:);
    end
end

end

%% ========================================================================
%  SINGLE PCC SECTION  (Frenet–Serret + geometric RK4)
%% ========================================================================

function [s, r, T, N, B] = rodSegment( ...
    curv_fun, tors_fun, Lseg, nPts, r0, T0, N0, B0)

s = linspace(0, Lseg, nPts)';

r = zeros(nPts, 3);
T = zeros(nPts, 3);
N = zeros(nPts, 3);
B = zeros(nPts, 3);

r(1,:) = r0; T(1,:) = T0; N(1,:) = N0; B(1,:) = B0;

for k = 1:nPts-1

    ds  = s(k+1) - s(k);
    kap = curv_fun(s(k));
    tau = tors_fun(s(k));

    y0 = [T(k,:)'; N(k,:)'; B(k,:)'];

    %% ---- Frenet–Serret right-hand side

    f = @(y) [ kap * y(4:6); ...
              -kap * y(1:3) + tau * y(7:9); ...
              -tau * y(4:6)];

    %% ---- RK4 steps

    k1 = f(y0);
    k2 = f(y0 + 0.5*ds*k1);
    k3 = f(y0 + 0.5*ds*k2);
    k4 = f(y0 +     ds*k3);

    y1 = y0 + (ds/6) * (k1 + 2*k2 + 2*k3 + k4);

    %% ---- re-orthonormalise frame

    T(k+1,:) = y1(1:3)' / norm(y1(1:3));
    N(k+1,:) = y1(4:6)' / norm(y1(4:6));
    B(k+1,:) = cross(T(k+1,:), N(k+1,:));
    B(k+1,:) = B(k+1,:) / norm(B(k+1,:));
    N(k+1,:) = cross(B(k+1,:), T(k+1,:));

    %% ---- integrate position

    r(k+1,:) = r(k,:) + ds * T(k,:);

end

end

%% ========================================================================
%  DRAW SPHERE  (solid, semi-transparent)
%% ========================================================================

function drawSphere(center, R)

[X, Y, Z] = sphere(60);

X = X * R + center(1);
Y = Y * R + center(2);
Z = Z * R + center(3);

surf(X, Y, Z, ...
    'FaceAlpha',  0.35, ...
    'EdgeColor',  'none', ...
    'FaceColor',  [0.9  0.3  0.2]);   % reddish obstacle

end

%% ========================================================================
%  DRAW SPHERE  (wireframe safety margin)
%% ========================================================================

function drawSphereWire(center, R)

[X, Y, Z] = sphere(30);

X = X * R + center(1);
Y = Y * R + center(2);
Z = Z * R + center(3);

surf(X, Y, Z, ...
    'FaceAlpha',  0.0, ...
    'EdgeColor',  [0.8  0.6  0.0], ...
    'EdgeAlpha',  0.25,            ...
    'LineStyle',  ':');

end