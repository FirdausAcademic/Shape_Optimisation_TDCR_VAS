%% Tip‐Position Control for a 3‐Section Continuum Robot

clc; clear; close all;

%% 1) Robot parameters
Ls    = [0.1,0.1,0.1];      % lengths of 3 sections [m]
L_total=sum(Ls);
nSec = numel(Ls);
r_d   = 0.009;  % disk radii
Ns    = 10*ones(1,nSec);         % disks per section


%% 2) Generalized Point Generation (with per-section limits)
trajectoryDensity = 4;
dt = 1/1000;

nSec = numel(Ls);

%── 2.1) Define per-section limits ───────────────────────────────────────
% Each row corresponds to one section: [minPull, maxPull]
lv_limits = [ ...
    0.001,  0.01;   % section 1
    0.001,  0.01; % section 2
    0.001,  0.01    % section 3
    ];

% Each row: [minAlpha, maxAlpha]
alpha_limits = [ ...
    0*pi,    0.5*pi;    % section 1
    0*pi,   1*pi;    % section 2
    0*pi,    2*pi    % section 3
    ];

%── 2.2) Normalized time vector ─────────────────────────────────────────
t = linspace(0, 1, trajectoryDensity);  % 1×N

%── 2.3) Build the shape-space trajectories ─────────────────────────────
% preallocate
del_lvstraj = zeros(nSec, trajectoryDensity);
alphastraj  = zeros(nSec, trajectoryDensity);

for i = 1:nSec
    % linear interpolation between [min, max] over t
    del_lvstraj(i,:) = linspace(lv_limits(i,1), lv_limits(i,2), trajectoryDensity);
    alphastraj(i,:)  = linspace(alpha_limits(i,1),alpha_limits(i,2),trajectoryDensity);
end

%%%% Preallocate
p_des = zeros(3, trajectoryDensity);

%%%% Build p_des in one pass
for k = 1:trajectoryDensity
    [x_d,~,~] = forwardKinematicsMultiSection(Ls, r_d,  del_lvstraj(:,k), alphastraj(:,k));
    p_des(:,k) = x_d;

end
%%
figure
plot3(p_des(1,:),p_des(2,:),p_des(3,:),'-g','LineWidth',2)
grid minor
span = 1.2*Ls(1);
xlim([-nSec*span, nSec*span]);  ylim([-nSec*span, nSec*span]);  zlim([-nSec*span, nSec*span]);

%% Optemization Start
% sections length 
Lstemp    = [0.1,0.1,0.1,0.1,0.1,0.1,0.1];


% Define a color palette (you can add more colors as needed)
colors = [
    1 0 0;    % red
    0 1 0;    % green
    0 0 1;    % blue
    1 0 1;    % magenta
    0 1 1;    % cyan
    1 1 0;    % yellow
    0.5 0 0;  % dark red
    0 0.5 0;  % dark green
    0 0 0.5;   % dark blue
    ];



for i=1:trajectoryDensity
    r_des=p_des(:,i);
    % for i=1
    %r_des=p_des0;
    figure
    % Initialize arrays to store plot handles and legend labels
    plotHandles = [];
    legendLabels = {};
    % Add their legend labels at the start
    legendLabels{1} = 'TDCR Tip';
    legendLabels{2} = 'Section Tip';
    h_des = plot3(r_des(2), r_des(3), r_des(1), '*r', 'MarkerSize', 8, 'LineWidth', 1.5);
    hold on;
    grid on;
    axis equal;
    box on
    xlabel('Y (m)', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel('Z (m)', 'Interpreter', 'latex', 'FontSize', 14);
    zlabel('X (m)', 'Interpreter', 'latex', 'FontSize', 14);
    span = 1.5 * Ls(1);
    xlim([-nSec*span, nSec*span]);
    ylim([-nSec*span, nSec*span]);
    zlim([-nSec*0, nSec*span]);
    for ii=1:length(Lstemp)-2
        % Get current color (cycling through palette)
        currentColor = colors(mod(ii, size(colors, 1)) + 1, :);

        % Create legend label (3 Sec, 4 Sec, etc.)
        legendLabels{end+1} = sprintf('%d Sections', 2+ii);

        % Plot rLine and store the handle
        Ls = Lstemp(1, 1:2+ii);


        m = numel(Ls);
        kappa0 = 2*ones(1,m);             % initial guess κ
        phi0   = zeros(1,m);              %  φ0 … φ2   (base + 2 junctions)
        lb_k = -6;  ub_k = 6;
        lb_p = -pi; ub_p =  pi;
        E_mod=210e9; I_area=1e-8; nu=0.3; G_mod=E_mod/(2*(1+nu)); J_area=2*I_area;

        [kOpt,phiOpt,Emin,err] = optimiseDescRod(Ls,r_des,kappa0,phi0, ...
            E_mod,I_area,G_mod,J_area,lb_k,ub_k,lb_p,ub_p,[]);

        fprintf('Tip err = %.2e m  |  E = %.4f J\n', err, Emin);
        [~,rLine,r_junc] = rodFK(kOpt,phiOpt,Ls,ceil(200*Ls/max(Ls)));
        %%




        % Store the plot handle


        h_junc = plot3(r_junc(:,2), r_junc(:,3), r_junc(:,1), '.g', 'MarkerSize', 20, 'LineWidth', 3);
        h_line = plot3(rLine(:,2), rLine(:,3), rLine(:,1), 'Color', currentColor, 'LineWidth', 2);
        plotHandles(end+1) = h_line;
        %plot3( [0, p_des0(2)], [0, p_des0(3)],[0, p_des0(1)], '--k', 'LineWidth', 1.5)
        legend([h_des, h_junc, plotHandles], legendLabels, 'Interpreter', 'latex', 'FontSize', 12);

        % Set the title
        %title('Minimum Energy Shape of TDCR', 'Interpreter', 'latex', 'FontSize', 14);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Add grid, axis labels, and limits




    end


end
% Create the legend with all handles and labels





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  DISCRETE-ROD TOOLBOX  (constant curvature + twist optimisation)
%  – with a BASE-TWIST φ0   (φ vector length = m)
% ------------------------------------------------------------------------
%  Save as e.g.  rod_cc_baseTwist.m   and run the demo at the bottom.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% =======================================================================
%  1)  MAIN OPTIMISER
% =======================================================================
function [kappa_opt, phi_opt, E_opt, tip_err] = optimiseDescRod( ...
    Ls, r_des, kappa0, phi0, ...                    % geometry
    E_mod, I_area, G_mod, J_area, ...              % material
    lb_kappa, ub_kappa, lb_phi, ub_phi, ...        % bounds
    nPts_vec)
% optimiseDescRod  Optimal {κ_i} & {φ_j} (φ₀ at base) to hit r_des.
%   φ vector length = m (φ₀ base  +  φ₁…φₘ₋₁ at junctions)
% -----------------------------------------------------------------------
%  OUTPUTS: kappa_opt, phi_opt, E_opt (J), tip_err (m)

% ---------- sizes & default sampling ----------
Ls = Ls(:)';   m = numel(Ls);
if nargin<13 || isempty(nPts_vec)
    nPts_vec = max(50, ceil(200*Ls/max(Ls)));
end

% ---------- variable pack ----------
x0 = [kappa0(:).'  ,  phi0(:).'];           % 1×(2m-1)
lb = [lb_kappa(:).', lb_phi(:).'];
ub = [ub_kappa(:).', ub_phi(:).'];

% ---------- objective ----------
obj = @(x) energyObjective(x, Ls, E_mod,I_area,G_mod,J_area);

% ---------- nonlinear equality (tip) ----------
nonlcon = @(x) tipConstraint(x, Ls, nPts_vec, r_des);

opts = optimoptions('fmincon','Algorithm','sqp', ...
    'Display','iter','MaxFunctionEvaluations',1e4);

[x_opt, E_opt] = fmincon(obj, x0, [],[],[],[], lb,ub, nonlcon, opts);

kappa_opt = x_opt(1:m);
phi_opt   = x_opt(m+1:end);

% tip error with final design
[~, r_line] = rodFK(kappa_opt, phi_opt, Ls, nPts_vec);
tip_err = norm(r_line(end,:) - r_des(:)');
end

%% =======================================================================
%  2)  Energy objective  (bending + all twists, incl. φ₀)
% =======================================================================
function E = energyObjective(x, Ls, E_mod,I_area,G_mod,J_area)
m = numel(Ls);
kappa = x(1:m);
phi   = x(m+1:end);        % length m

E_bend  = 0.5*E_mod*I_area * sum( Ls .* kappa.^2 );
E_twist = 0.5*G_mod*J_area * sum( phi.^2 );
E = E_bend + E_twist;
end

%% =======================================================================
%  3)  Non-linear equality: tip coincidence
% =======================================================================
function [c,ceq] = tipConstraint(x, Ls, nPts_vec, r_des)
m = numel(Ls);
kappa = x(1:m);
phi   = x(m+1:end);
[~, r_line] = rodFK(kappa, phi, Ls, nPts_vec);
ceq = r_line(end,:).' - r_des(:);
c   = [];
end

%% =======================================================================

function [s_tot, r_tot, r_junc] = rodFK(kappa, phi, Ls, nPts_vec)
% rodFK  Forward kinematics for an m-segment constant-curvature rod
%        with base twist φ0 and junction twists φ1…φ_{m−1}.
%
% Inputs
%   kappa     : 1×m  constant curvature of each segment       [1/m]
%   phi       : 1×m  twists  [φ0 (base), φ1 … φ_{m−1}]        [rad]
%   Ls        : 1×m  lengths of segments                      [m]
%   nPts_vec  : 1×m  sample points per segment   (≥2 each)
%
% Outputs
%   s_tot     : column vector of arc-length samples            [m]
%   r_tot     : N×3  centre-line coordinates                   [m]
%   r_junc    : (m+1)×3 coordinates of base, all junctions, tip[m]

% ------------------------------------------------------------
% 0)  basic sizes & helpers
% ------------------------------------------------------------
m = numel(kappa);
if numel(phi) ~= m
    error('phi must have length m  (φ0 base  +  φ1..φ_{m-1})');
end
if numel(Ls) ~= m || numel(nPts_vec) ~= m
    error('Ls and nPts_vec must also have length m');
end

% ------------------------------------------------------------
% 1)  build constant-strain function handles
% ------------------------------------------------------------
curvature_funs = cellfun(@(k)@(s)k, num2cell(kappa),  'uni',0);
torsion_funs   = repmat({@(s)0}, 1, m);  % still planar per segment

% ------------------------------------------------------------
% 2)  initial frame  (apply base twist φ0 to N,B)
% ------------------------------------------------------------
T0   = [1 0 0];                 % tangent at the base
Nor0 = [0 1 0]*cos(phi(1)) + [0 0 1]*sin(phi(1));
B0   = cross(T0, Nor0);

% ------------------------------------------------------------
% 3)  run multi-section integrator (junction twists = φ1..φ_{m-1})
% ------------------------------------------------------------
[s_tot, r_tot, ~, ~, ~] = rod_multi_sections_customStart( ...
    curvature_funs, torsion_funs, ...
    Ls, nPts_vec,                                   ...
    phi(2:end),                                     ... % junction twists
    [0 0 0], T0, Nor0, B0);                         % custom start frame

% ------------------------------------------------------------
% 4)  extract coordinates of every junction
% ------------------------------------------------------------
cumL = [0  cumsum(Ls)];          % 0, L1, L1+L2, … Ltot   (length m+1)
r_junc = zeros(numel(cumL), 3);

for j = 1:numel(cumL)
    [~, idx] = min( abs( s_tot - cumL(j) ) );  % nearest sampled point
    r_junc(j,:) = r_tot(idx,:);
end
end

%% =======================================================================
%  5)  Multi-section integrator with custom start frame
% =======================================================================
function [s_tot,r_tot,T_tot,Nor_tot,B_tot] = rod_multi_sections_customStart( ...
    curvature_funs,torsion_funs,Ls,nPts_vec,phis, ...
    r0,T0,N0,B0)
% identical to your old rod_multi_sections but uses (r0,T0,N0,B0).
m = numel(Ls);    s_tot=[]; r_tot=[]; T_tot=[]; Nor_tot=[]; B_tot=[];
r_prev=r0; T_prev=T0; Nor_prev=N0; B_prev=B0; cumL=0;
for i=1:m
    Li=Ls(i); ni=nPts_vec(i);
    [s_i,r_i,T_i,Nor_i,B_i]=rod_segment( ...
        curvature_funs{i},torsion_funs{i},Li,ni,r_prev,T_prev,Nor_prev,B_prev);
    s_i_shift=cumL+s_i;
    if i==1
        s_tot=s_i_shift; r_tot=r_i; T_tot=T_i; Nor_tot=Nor_i; B_tot=B_i;
    else
        s_tot=[s_tot; s_i_shift(2:end)];
        r_tot=[r_tot; r_i(2:end,:)]; T_tot=[T_tot; T_i(2:end,:)];
        Nor_tot=[Nor_tot; Nor_i(2:end,:)]; B_tot=[B_tot; B_i(2:end,:)];
    end
    cumL=cumL+Li;
    if i<m
        phi=phis(i);
        T_prev=T_i(end,:); Nend=Nor_i(end,:); Bend=B_i(end,:);
        Nor_prev= cos(phi)*Nend+sin(phi)*Bend;
        B_prev  =-sin(phi)*Nend+cos(phi)*Bend;
        r_prev=r_i(end,:);
    end
end
end

%% =======================================================================
%  6)  Single-segment RK4 integrator (unchanged)
% =======================================================================
function [s_vec,r,T,Nor,B] = rod_segment(curv_fun,tau_fun,L_seg,nPts,r0,T0,N0,B0)
s_vec=linspace(0,L_seg,nPts).';
r=zeros(nPts,3);
r(1,:)=r0;
T=zeros(nPts,3);
T(1,:)=T0;
Nor=zeros(nPts,3);
Nor(1,:)=N0;
B=zeros(nPts,3);
B(1,:)=B0;
for k=1:nPts-1
    ds=s_vec(k+1)-s_vec(k); s=s_vec(k);
    kappa=curv_fun(s); tau=tau_fun(s);
    y0=[T(k,:)'; Nor(k,:)'; B(k,:)'];
    f=@(y)[kappa*y(4:6); -kappa*y(1:3)+tau*y(7:9); -tau*y(4:6)];
    k1=f(y0);
    k2=f(y0+0.5*ds*k1);
    k3=f(y0+0.5*ds*k2);
    k4=f(y0+ds*k3);
    y1=y0+(ds/6)*(k1+2*k2+2*k3+k4);
    T(k+1,:)=y1(1:3)'/norm(y1(1:3));
    Nor(k+1,:)=y1(4:6)'/norm(y1(4:6));
    B(k+1,:)=cross(T(k+1,:),Nor(k+1,:));
    B(k+1,:)=B(k+1,:)/norm(B(k+1,:));
    Nor(k+1,:)=cross(B(k+1,:),T(k+1,:));
    r(k+1,:)=r(k,:)+ds*T(k,:);
end
end
