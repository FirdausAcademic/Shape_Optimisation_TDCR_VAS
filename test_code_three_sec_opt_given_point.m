clear;  clc;
 close all;

% ----- rod layout -----
Ls   = [0.3 0.25 0.35];      m = numel(Ls);

% ----- material -----
E_mod  = 210e9;   I_area = 1e-8;
nu     = 0.3;     G_mod  = E_mod/(2*(1+nu));
J_area = 2*I_area;

% ----- target tip position -----
% r_des = [0.2649 0.5251 0.3627];
% r_des = [0.1732 0.5979 0.2456];
% r_des = [0.3835 0.6458 0.2278];
r_des = [0.7403 -0.3626 0.2278];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Optemisation
% ----- optimiser setup -----
kappa0 = 2*ones(1,m);          % initial κ
phi0   = zeros(1,m-1);         % initial φ
lb_k   = -6*ones(1,m);  ub_k = 6*ones(1,m);
lb_p   = -pi*ones(1,m-1); ub_p = pi*ones(1,m-1);

% (optional) sampling per segment
nPts_vec = [];                 % [] → use default inside the function

% ----- call optimiser -----
[kappa_opt, phi_opt, E_opt] = optimize_desc_rod( ...
        Ls, r_des, kappa0, phi0, ...
        E_mod, I_area, G_mod, J_area, ...
        lb_k, ub_k, lb_p, ub_p, ...
        nPts_vec );

% ----- report -----
fprintf('\\nOptimum found\\n');
disp(table((1:m)', kappa_opt', 'VariableNames',{'seg','kappa'}));
disp(table((1:m-1)', phi_opt', 'VariableNames',{'junction','phi'}));
fprintf('Energy     : %.4f J\\n', E_opt);
% fprintf('Tip error  : %.2e m\\n',error );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% kinematics
% define 3 segments

curvature_funs = { @(s)kappa_opt(1), @(s) kappa_opt(2), @(s)kappa_opt(3)};
torsion_funs   = { @(s)0, @(s)0, @(s)0 };           % planar
nPts_vec       = [200,        150,           250];
phis           = [phi_opt(1),      phi_opt(2)];               % two junctions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 2) Material / cross-section constants
% E_mod  = 210e9;                % Young's modulus  [Pa]
% I_area = 1.0e-8;               % 2nd moment       [m^4]
% nu     = 0.3;                  % Poisson ratio
% G_mod  = E_mod/(2*(1+nu));     % shear modulus    [Pa]
% J_area = 2*I_area;             % polar moment     [m^4] (circular)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[s_tot, r_tot, ~,~,~] = rod_multi_sections( ...
    curvature_funs, torsion_funs, Ls, nPts_vec, phis );

figure;
plot3(r_tot(:,1), r_tot(:,2), r_tot(:,3),'LineWidth',1.6);
hold on
plot3(r_des(1),r_des(2),r_des(3),'*r','LineWidth',3,'MarkerSize',6)
axis equal; grid on;
%plot3(r_des(1),r_des(2),r_des(3),'*r','LineWidth',4,'MarkerSize',8)
xlabel('x'); ylabel('y'); zlabel('z');
title('m-segment rod with discrete twists at the junctions');


%% Function
function [kappa_opt, phi_opt, E_opt] = optimize_desc_rod( ...
            Ls, r_des, kappa0, phi0, ...
            E_mod, I_area, G_mod, J_area, ...
            lb_kappa, ub_kappa, lb_phi, ub_phi, ...
            nPts_vec)
% OPTIMIZE_PLANAR_CC  constant-curvature + discrete-twist optimiser
%
%   [kappa_opt, phi_opt, E_opt, tip_err] = optimize_planar_cc(...)
%
% Inputs (vectors are row or column OK)
%   Ls       : 1×m   segment lengths
%   r_des    : 3×1   desired tip position
%   kappa0   : 1×m   initial guess for κ
%   phi0     : 1×(m-1) initial guess for φ
%   E_mod,I_area,G_mod,J_area : material constants
%   lb_kappa,ub_kappa         : bounds on κ  (scalar or 1×m)
%   lb_phi,ub_phi             : bounds on φ  (scalar or 1×(m-1))
%   nPts_vec : sampling per segment (default ceil(200*Ls/max(Ls)))
%
% Outputs
%   kappa_opt : optimal curvatures [1×m]
%   phi_opt   : optimal twists     [1×(m-1)]
%   E_opt     : minimum energy     [J]
%   tip_err   : ‖r_tip - r_des‖    [m]

Ls       = Ls(:)';
m        = numel(Ls);
if nargin<13 || isempty(nPts_vec)     % << changed 12 → 13
    nPts_vec = ceil(200*Ls/max(Ls));
end
    %% pack initial vector  x = [κ1..κm  φ1..φ_{m-1}]
    x0  = [kappa0(:)'  ,  phi0(:)'];
    lb  = [lb_kappa(:)' , lb_phi(:)'];
    ub  = [ub_kappa(:)' , ub_phi(:)'];

    %% objective (energy) ----------------------------------------------
    function E = objFun(x)
        kappa = x(1:m);
        phi   = x(m+1:end);
        % build constant curvature/torsion handles
        curvature_funs = cellfun(@(k)@(s)k, num2cell(kappa),'uni',0);
        torsion_funs   = repmat({@(s)0}, 1, m);      % planar
        E = rod_multi_sections_energy(curvature_funs, torsion_funs, ...
                    Ls, nPts_vec, phi, E_mod, I_area, G_mod, J_area);
    end

%% nonlinear equality: tip-pos constraint --------------------------
function [c,ceq] = nonlcon(x)
    kappa = x(1:m);
    phi   = x(m+1:end);

    % build constant-strain function handles
    curvature_funs = cellfun(@(k)@(s)k, num2cell(kappa),'uni',0);
    torsion_funs   = repmat({@(s)0}, 1, m);   % zero torsion per segment

    % full 3-D forward kinematics (same one you plot with)
    [~, r_tot] = rod_multi_sections(curvature_funs, torsion_funs, ...
                                    Ls, nPts_vec, phi);

    r_tip = r_tot(end,:);

    c   = [];                        % no inequality constraints
    ceq = r_tip(:) - r_des(:);       % force x,y,z to match target
end

    opts = optimoptions('fmincon', ...
            'Algorithm','sqp', ...
            'Display','iter', ...
            'MaxFunctionEvaluations',1e4);

    [x_opt, E_opt] = fmincon(@objFun, x0, [],[],[],[], lb,ub, @nonlcon, opts);

    kappa_opt = x_opt(1:m);
    phi_opt   = x_opt(m+1:end);

end


function [E_total, E_bend, E_tors, E_twist] = rod_multi_sections_energy( ...
        curvature_funs, torsion_funs, Ls, nPts_vec, phis, ...
        E_mod, I_area, G_mod, J_area)
% ROD_MULTI_SECTIONS_ENERGY  Total elastic energy for an m‑segment rod
% with optional discrete twist (phi) at each junction.
%
%   [E_total, E_bend, E_tors, E_twist] = rod_multi_sections_energy(...)
%
% Inputs
%   curvature_funs : 1×m cell array of @(s)->kappa_i(s)
%   torsion_funs   : 1×m cell array of @(s)->tau_i(s)
%   Ls             : 1×m vector of segment lengths [m]
%   nPts_vec       : 1×m vector of sample counts per segment
%   phis           : 1×(m‑1) vector of twist angles (rad) at the m‑1 junctions.
%                     Use zeros(m-1,1) if no twist.  Positive φ follows right‑hand rule
%                     about the downstream tangent.
%   E_mod, I_area  : Young's modulus & second moment of area
%   G_mod, J_area  : shear modulus & polar moment of area (for circular / isotropic rod)
%
% Outputs
%   E_total  : scalar, total elastic energy  (bending + torsion + twist) [J]
%   E_bend   : 1×m vector, bending energy of each segment               [J]
%   E_tors   : 1×m vector, torsional energy of each segment             [J]
%   E_twist  : 1×(m‑1) vector, discrete twist energy at each junction   [J]
%
% Notes
%   • Bending & torsional energies follow the usual distributed formulae
%       ½∫ EI κ² ds   and   ½∫ GJ τ² ds
%   • Discrete twist energy is modelled as a torsional spring located at
%     the junction:   U_i = ½ GJ φ_i²   (dimensionally identical to
%     the distributed expression if we imagine the twist localised over a
%     vanishingly small length).
%   • All variable names kept consistent with previous scripts.

    m = numel(Ls);
    assert(numel(curvature_funs)==m && numel(torsion_funs)==m, "Mismatch in segment counts");
    assert(numel(nPts_vec)==m, "nPts_vec length must equal number of segments");
    if isempty(phis); phis = zeros(1,m-1); end
    assert(numel(phis)==m-1, "phis must have m‑1 elements");

    % Pre‑allocate energy arrays
    E_bend  = zeros(1,m);
    E_tors  = zeros(1,m);
    E_twist = 0.5 * G_mod * J_area * phis.^2;   % junction energies

    % Loop over segments and integrate distributed energies
    for i = 1:m
        Li = Ls(i);
        ni = nPts_vec(i);
        [~, E_bend(i), E_tors(i)] = computeRodEnergy( ...
            curvature_funs{i}, torsion_funs{i}, Li, ni, ...
            E_mod, I_area, G_mod, J_area );
    end

    % Total
    E_total = sum(E_bend) + sum(E_tors) + sum(E_twist);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Helper : computeRodEnergy (unchanged from earlier scripts)            %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [E_total, E_bend, E_tors] = computeRodEnergy(curvature_fun, torsion_fun, L, nPts, E_mod, I_area, G_mod, J_area)
    s     = linspace(0, L, nPts).';              % column
    kappa = curvature_fun(s);   tau = torsion_fun(s);
    if isscalar(kappa), kappa = kappa * ones(size(s)); else, kappa = kappa(:); end
    if isscalar(tau),   tau   = tau   * ones(size(s)); else, tau   = tau(:);   end
    E_bend = 0.5 * E_mod * I_area * trapz(s, kappa.^2);
    E_tors = 0.5 * G_mod * J_area * trapz(s, tau  .^2);
    E_total = E_bend + E_tors;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% forward Kinematics 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [s_tot, r_tot, T_tot, Nor_tot, B_tot] = rod_multi_sections( ...
        curvature_funs, torsion_funs, Ls, nPts_vec, phis)
% ROD_MULTI_SECTIONS  m‐segment rod with discrete twist at each junction
%
% Inputs
%   curvature_funs : 1×m cell of @(s)->κ_i(s)
%   torsion_funs   : 1×m cell of @(s)->τ_i(s)
%   Ls             : 1×m vector of segment lengths [L1,…,Lm]
%   nPts_vec       : 1×m vector of sample counts per segment
%   phis           : 1×(m-1) vector of twist angles [φ1,…,φ_{m−1}]
%
% Outputs (all concatenated length = sum(nPts_vec)−(m−1))
%   s_tot   : cumulative arc‐length
%   r_tot   : centre‐line positions
%   T_tot   : tangents
%   Nor_tot : normals
%   B_tot   : binormals

    m = numel(Ls);
    assert( numel(curvature_funs)==m && numel(torsion_funs)==m );
    assert( numel(nPts_vec)==m && numel(phis)==m-1 );
    
    % initial frame at base
    r_prev   = [0,0,0];
    T_prev   = [1,0,0];
    Nor_prev = [0,1,0];
    B_prev   = cross(T_prev, Nor_prev);
    
    s_tot = []; r_tot = []; T_tot = []; Nor_tot = []; B_tot = [];
    cumL = 0;
    
    for i = 1:m
        Li  = Ls(i);
        ni  = nPts_vec(i);
        % integrate section i
        [s_i, r_i, T_i, Nor_i, B_i] = rod_segment( ...
            curvature_funs{i}, torsion_funs{i}, ...
            Li, ni, r_prev, T_prev, Nor_prev, B_prev );
        
        % shift s_i and stitch (drop duplicate start for i>1)
        s_i_shift = cumL + s_i;
        if i==1
            s_tot   = s_i_shift;
            r_tot   = r_i;
            T_tot   = T_i;
            Nor_tot = Nor_i;
            B_tot   = B_i;
        else
            s_tot   = [s_tot;   s_i_shift(2:end)];
            r_tot   = [r_tot;   r_i(2:end,:)];
            T_tot   = [T_tot;   T_i(2:end,:)];
            Nor_tot = [Nor_tot; Nor_i(2:end,:)];
            B_tot   = [B_tot;   B_i(2:end,:)];
        end
        
        % prepare for next segment
        cumL = cumL + Li;
        if i < m
            % discrete twist φ_i about T_i(end)
            phi   = phis(i);
            T_prev   = T_i(end,:);
            Nend     = Nor_i(end,:);
            Bend     = B_i(end,:);
            Nor_prev =  cos(phi)*Nend + sin(phi)*Bend;
            B_prev   = -sin(phi)*Nend + cos(phi)*Bend;
            r_prev   = r_i(end,:);
        end
    end
end
function [s_vec, r, T, Nor, B] = rod_segment( ...
        curvature_fun, torsion_fun, L_seg, nPts, r0, T0, Nor0, B0)

    s_vec = linspace(0, L_seg, nPts).';
    r   = zeros(nPts,3);  r(1,:)   = r0;
    T   = zeros(nPts,3);  T(1,:)   = T0;
    Nor = zeros(nPts,3);  Nor(1,:) = Nor0;
    B   = zeros(nPts,3);  B(1,:)   = B0;

    for k = 1:nPts-1
        ds = s_vec(k+1) - s_vec(k);
        s  = s_vec(k);
        % get local curvature & torsion
        kappa = curvature_fun(s);
        tau   = torsion_fun(s);

        y0 = [T(k,:)'; Nor(k,:)'; B(k,:)'];
        f  = @(y)[ ...
            kappa * y(4:6);                % dT/ds   = κ·Nor
           -kappa * y(1:3) + tau * y(7:9);  % dNor/ds = -κ·T + τ·B
           -tau   * y(4:6) ];              % dB/ds   = -τ·Nor

        % RK4
        k1 = f(y0);
        k2 = f(y0 + 0.5*ds*k1);
        k3 = f(y0 + 0.5*ds*k2);
        k4 = f(y0 +    ds*k3);
        y1 = y0 + (ds/6)*(k1 + 2*k2 + 2*k3 + k4);

        % re‐orthonormalise
        T(k+1,:)   = y1(1:3)'   / norm(y1(1:3));
        Nor(k+1,:) = y1(4:6)'   / norm(y1(4:6));
        B(k+1,:)   = cross(T(k+1,:), Nor(k+1,:));
        B(k+1,:)   = B(k+1,:)   / norm(B(k+1,:));
        Nor(k+1,:) = cross(B(k+1,:), T(k+1,:));

        r(k+1,:) = r(k,:) + ds * T(k,:);
    end
end
