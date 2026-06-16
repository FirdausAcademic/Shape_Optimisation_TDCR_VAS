clc;
clear;
close all;

%% ========================================================================
%  TDCR PCC OPTIMISATION
%  ---------------------------------------------------------------
%  METHOD:
%
%   TDCR Piecewise Constant Curvature
% + Frenet Frame Integration
% + Geometric RK4 Integration
% + Log Barrier Cylinder Avoidance
% + Augmented Lagrangian Tip Targeting
%
%  Unknowns:
%       x = [kappa_1 ... kappa_m  phi_1 ... phi_m]
%
%  Optimisation:
%
%     F(x) =
%       E_elastic
%     + E_wall
%     + E_tip(AL)
%
%  Uses:
%       fminunc (unconstrained optimisation)
%
%% ========================================================================

%% ========================================================================
%  TDCR PARAMETERS
%% ========================================================================

Ls = [0.10 0.10 0.10 0.10 0.10];

m = numel(Ls);

Ltot = sum(Ls);

%% ========================================================================
%  TARGET TIP
%% ========================================================================

r_des = [0.07 0.04 0.38];

%% ========================================================================
%  MATERIAL
%% ========================================================================

EI = 1.0;
GJ = 0.1;

%% ========================================================================
%  CYLINDER
%% ========================================================================

cyl.radius = 0.06;
cyl.length = 0.25;

cyl.axis_point = [0 0 0];

cyl.axis_dir = [0 0 1];
cyl.axis_dir = cyl.axis_dir / norm(cyl.axis_dir);

d_safe = 0.003;

%% ========================================================================
%  BARRIER + AL PARAMETERS
%% ========================================================================

mu_wall   = 1e-1;

mu_tip    = 10;

mu_scale  = 5;

lambda    = zeros(3,1);

nAL       = 8;

%% ========================================================================
%  DISCRETISATION
%% ========================================================================

nPts = 80 * ones(1,m);

%% ========================================================================
%  INITIAL GUESS
%% ========================================================================

x0 = [2*ones(1,m)  zeros(1,m)];

x = x0;

%% ========================================================================
%  OPTIMISER
%% ========================================================================

opts = optimoptions('fminunc',...
    'Algorithm','quasi-newton',...
    'Display','iter',...
    'MaxIterations',500,...
    'MaxFunctionEvaluations',2e5,...
    'OptimalityTolerance',1e-8,...
    'StepTolerance',1e-10);

%% ========================================================================
%  AUGMENTED LAGRANGIAN OUTER LOOP
%% ========================================================================

fprintf('\n');
fprintf('=====================================================\n');
fprintf(' TDCR PCC + Barrier + Augmented Lagrangian\n');
fprintf('=====================================================\n\n');

for al = 1:nAL

    fprintf('AL ITERATION %d\n',al);

    obj = @(xx) totalCost(...
        xx,...
        Ls,...
        nPts,...
        EI,...
        GJ,...
        cyl,...
        d_safe,...
        r_des,...
        mu_wall,...
        lambda,...
        mu_tip);

    x = fminunc(obj,x,opts);

    %% ---- evaluate tip error

    [~,rline] = rodFK(...
        x(1:m),...
        x(m+1:end),...
        Ls,...
        nPts);

    e = rline(end,:)' - r_des(:);

    fprintf('Tip Error = %.6e\n',norm(e));

    %% ---- update AL

    lambda = lambda + mu_tip * e;

    mu_tip = mu_tip * mu_scale;

    fprintf('mu_tip    = %.3e\n',mu_tip);

    fprintf('--------------------------------------------\n');

    if norm(e) < 1e-5
        break;
    end
end

%% ========================================================================
%  FINAL RESULT
%% ========================================================================

kappa = x(1:m);

phi = x(m+1:end);

[s,rline,rj] = rodFK(kappa,phi,Ls,nPts);

Eel = elasticEnergy(kappa,phi,Ls,EI,GJ);

fprintf('\n');
fprintf('=====================================================\n');
fprintf('FINAL RESULT\n');
fprintf('=====================================================\n');

fprintf('Elastic Energy = %.6f\n',Eel);

fprintf('Tip Error       = %.6e\n',...
    norm(rline(end,:)-r_des));

fprintf('\n');

%% ========================================================================
%  VISUALISATION
%% ========================================================================

figure('Color','w');

hold on;
grid on;
axis equal;
view(3);

xlabel('X');
ylabel('Y');
zlabel('Z');

title('TDCR PCC Optimisation Inside Cylinder');

%% ---- cylinder

drawCylinder(cyl.radius,cyl.length);

%% ---- rod

plot3(...
    rline(:,1),...
    rline(:,2),...
    rline(:,3),...
    'r-','LineWidth',3);

%% ---- junctions

plot3(...
    rj(:,1),...
    rj(:,2),...
    rj(:,3),...
    '.k','MarkerSize',25);

%% ---- target

plot3(...
    r_des(1),...
    r_des(2),...
    r_des(3),...
    'bp',...
    'MarkerSize',16,...
    'MarkerFaceColor','b');

legend(...
    'Cylinder',...
    'Rod',...
    'Junctions',...
    'Target');

%% ========================================================================
%  TOTAL COST
%% ========================================================================

function F = totalCost(...
    x,...
    Ls,...
    nPts,...
    EI,...
    GJ,...
    cyl,...
    d_safe,...
    r_des,...
    mu_wall,...
    lambda,...
    mu_tip)

m = numel(Ls);

kappa = x(1:m);

phi = x(m+1:end);

%% ---- elastic energy

Eel = elasticEnergy(kappa,phi,Ls,EI,GJ);

%% ---- forward kinematics

[~,rline] = rodFK(kappa,phi,Ls,nPts);

%% ---- wall barrier

Ewall = wallBarrier(...
    rline,...
    cyl,...
    d_safe,...
    mu_wall);

%% ---- tip AL

e = rline(end,:)' - r_des(:);

Etip = lambda' * e ...
    + 0.5 * mu_tip * (e' * e);

%% ---- total

F = Eel + Ewall + Etip;

fprintf(...
    'Eel = %.4e   Ewall = %.4e   Etip = %.4e\r',...
    Eel,Ewall,Etip);

end

%% ========================================================================
%  ELASTIC ENERGY
%% ========================================================================

function E = elasticEnergy(kappa,phi,Ls,EI,GJ)

Ebend = 0.5 * EI * sum(Ls .* kappa.^2);

Etwist = 0.5 * GJ * sum(phi.^2);

E = Ebend + Etwist;

end

%% ========================================================================
%  LOG BARRIER CYLINDER
%% ========================================================================

function Ewall = wallBarrier(...
    rline,...
    cyl,...
    d_safe,...
    mu_wall)

N = size(rline,1);

a = cyl.axis_dir(:)';

c0 = cyl.axis_point(:)';

R = cyl.radius;

Ewall = 0;

for i = 1:N

    p = rline(i,:);

    v = p - c0;

    t = dot(v,a);

    %% ---- only inside cylinder length

    if t >= 0 && t <= cyl.length

        vr = v - t*a;

        rho = norm(vr);

        gap = R - rho - d_safe;

        %% ---- collision

        if gap <= 0

            Ewall = Ewall + 1e6 * abs(gap)^2;

        else

            Ewall = Ewall ...
                - mu_wall * log(gap / R);

        end
    end
end

end

%% ========================================================================
%  FORWARD KINEMATICS
%% ========================================================================

function [s_tot,r_tot,r_junc] = rodFK(...
    kappa,...
    phi,...
    Ls,...
    nPts)

m = numel(kappa);

curv_funs = cellfun(...
    @(k) @(s) k,...
    num2cell(kappa),...
    'UniformOutput',false);

tors_funs = repmat({@(s)0},1,m);

%% ---- initial frame

T0 = [0 0 1];

N0 = [1 0 0]*cos(phi(1)) ...
    + [0 1 0]*sin(phi(1));

B0 = cross(T0,N0);

[s_tot,r_tot,~,~,~] = rodMulti(...
    curv_funs,...
    tors_funs,...
    Ls,...
    nPts,...
    phi(2:end),...
    [0 0 0],...
    T0,...
    N0,...
    B0);

%% ---- junctions

cumL = [0 cumsum(Ls)];

r_junc = zeros(numel(cumL),3);

for j = 1:numel(cumL)

    [~,idx] = min(abs(s_tot-cumL(j)));

    r_junc(j,:) = r_tot(idx,:);

end

end

%% ========================================================================
%  MULTI-SECTION PROPAGATION
%% ========================================================================

function [s_tot,r_tot,T_tot,N_tot,B_tot] = ...
    rodMulti(...
    curv_funs,...
    tors_funs,...
    Ls,...
    nPts,...
    phis,...
    r0,...
    T0,...
    N0,...
    B0)

m = numel(Ls);

s_tot=[];
r_tot=[];
T_tot=[];
N_tot=[];
B_tot=[];

r_prev=r0;
T_prev=T0;
N_prev=N0;
B_prev=B0;

cumL=0;

for i=1:m

    Li=Ls(i);

    ni=nPts(i);

    [s_i,r_i,T_i,N_i,B_i] = ...
        rodSegment(...
        curv_funs{i},...
        tors_funs{i},...
        Li,...
        ni,...
        r_prev,...
        T_prev,...
        N_prev,...
        B_prev);

    s_i = s_i + cumL;

    if i==1

        s_tot=s_i;
        r_tot=r_i;
        T_tot=T_i;
        N_tot=N_i;
        B_tot=B_i;

    else

        s_tot=[s_tot;s_i(2:end)];

        r_tot=[r_tot;r_i(2:end,:)];

        T_tot=[T_tot;T_i(2:end,:)];

        N_tot=[N_tot;N_i(2:end,:)];

        B_tot=[B_tot;B_i(2:end,:)];

    end

    cumL = cumL + Li;

    %% ---- rotate bending plane

    if i<m

        ph = phis(i);

        T_prev = T_i(end,:);

        Nend = N_i(end,:);

        Bend = B_i(end,:);

        N_prev = ...
            cos(ph)*Nend ...
            + sin(ph)*Bend;

        B_prev = ...
            -sin(ph)*Nend ...
            + cos(ph)*Bend;

        r_prev = r_i(end,:);

    end
end

end

%% ========================================================================
%  SINGLE PCC SECTION
%% ========================================================================

function [s,r,T,N,B] = ...
    rodSegment(...
    curv_fun,...
    tors_fun,...
    Lseg,...
    nPts,...
    r0,...
    T0,...
    N0,...
    B0)

s = linspace(0,Lseg,nPts)';

r=zeros(nPts,3);

T=zeros(nPts,3);

N=zeros(nPts,3);

B=zeros(nPts,3);

r(1,:)=r0;

T(1,:)=T0;

N(1,:)=N0;

B(1,:)=B0;

for k=1:nPts-1

    ds = s(k+1)-s(k);

    kap = curv_fun(s(k));

    tau = tors_fun(s(k));

    y0 = [...
        T(k,:)';...
        N(k,:)';...
        B(k,:)'];

    %% ---- Frenet equations

    f = @(y)[...
        kap*y(4:6);...
       -kap*y(1:3)+tau*y(7:9);...
       -tau*y(4:6)];

    %% ---- RK4

    k1=f(y0);

    k2=f(y0+0.5*ds*k1);

    k3=f(y0+0.5*ds*k2);

    k4=f(y0+ds*k3);

    y1 = y0 + (ds/6)*(...
        k1 + 2*k2 + 2*k3 + k4);

    %% ---- orthonormalise

    T(k+1,:) = y1(1:3)' ...
        / norm(y1(1:3));

    N(k+1,:) = y1(4:6)' ...
        / norm(y1(4:6));

    B(k+1,:) = cross(...
        T(k+1,:),...
        N(k+1,:));

    B(k+1,:) = ...
        B(k+1,:) ...
        / norm(B(k+1,:));

    N(k+1,:) = cross(...
        B(k+1,:),...
        T(k+1,:));

    %% ---- integrate position

    r(k+1,:) = ...
        r(k,:) + ds*T(k,:);

end

end

%% ========================================================================
%  DRAW CYLINDER
%% ========================================================================

function drawCylinder(R,H)

[X,Y,Z] = cylinder(R,80);

Z = Z * H;

surf(...
    X,Y,Z,...
    'FaceAlpha',0.5,...
    'EdgeColor','none',...
    'FaceColor',[0.7 0.8 1]);

end