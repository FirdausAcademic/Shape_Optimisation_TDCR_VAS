% % % clc; clear; close all;
% % % 
% % % % =========================================================================
% % % %  TDCR SHAPE OPTIMISATION — LOG-BARRIER OBSTACLE AVOIDANCE
% % % %  -----------------------------------------------------------------------
% % % %  This script merges the two previous approaches:
% % % %
% % % %   - Obstacle SET (geometry, plotting) is taken from the "general surface
% % % %     avoidance" script: sphere / plane / cylinder / ellipsoid obstacles,
% % % %     each described by a signed-distance function and a required
% % % %     clearance d_safe.
% % % %
% % % %   - Optimisation THEORY is taken from the "log-barrier" script:
% % % %       F(x) = E_elastic  +  sum_p E_barrier(obstacle p)  +  E_tip(AL)
% % % %     solved with an unconstrained solver (fminunc) and an Augmented
% % % %     Lagrangian (AL) outer loop that drives the tip to r_des while the
% % % %     log barrier keeps the backbone outside every obstacle's safety
% % % %     shell. The barrier blows up smoothly as the rod approaches any
% % % %     obstacle surface, giving well-behaved gradients (no hard min()).
% % % %
% % % %  Unknowns:  x = [kappa_1 ... kappa_m   phi_1 ... phi_m]
% % % %
% % % %  Each obstacle is a struct (same convention as the surface-avoidance
% % % %  script):
% % % %    .type     – 'plane' | 'sphere' | 'cylinder' | 'ellipsoid'
% % % %    .d_safe   – required clearance [m]
% % % %    .normal   – unit normal (plane)
% % % %    .point    – point on plane / centre of sphere,cylinder,ellipsoid
% % % %    .radius   – radius (sphere, cylinder)
% % % %    .axis     – unit axis direction (cylinder)
% % % %    .semiaxes – [a b c] (ellipsoid)
% % % %    .axes     – 3x3 orientation matrix, rows = principal directions (ellipsoid)
% % % % =========================================================================
% % % 
% % % %% -----------------------------------------------------------------------
% % % %  ROD PARAMETERS
% % % % -----------------------------------------------------------------------
% % % Ls    = [0.1, 0.15, 0.1, 0.15, 0.1];   % section lengths [m]
% % % m     = numel(Ls);
% % % L_tot = sum(Ls);
% % % 
% % % %% -----------------------------------------------------------------------
% % % %  TARGET TIP POSITION
% % % % -----------------------------------------------------------------------
% % % r_des = [0.2  0.15  0.15];
% % % 
% % % fprintf('Total rod length   = %.3f m\n', L_tot);
% % % fprintf('||r_des||          = %.3f m\n', norm(r_des));
% % % fprintf('Reachable?         = %s\n\n', mat2str(norm(r_des) <= L_tot));
% % % 
% % % %% -----------------------------------------------------------------------
% % % %  MATERIAL & CROSS-SECTION
% % % % -----------------------------------------------------------------------
% % % E_mod  = 210e9;
% % % I_area = 1e-8;
% % % nu     = 0.3;
% % % G_mod  = E_mod / (2*(1+nu));
% % % J_area = 2 * I_area;
% % % 
% % % EI = E_mod * I_area;     % bending stiffness  [N.m^2]
% % % GJ = G_mod * J_area;     % torsional stiffness [N.m^2]
% % % 
% % % %% -----------------------------------------------------------------------
% % % %  OBSTACLE DEFINITIONS  (same physical set as the surface-avoidance script)
% % % % -----------------------------------------------------------------------
% % % obs = {};
% % % 
% % % % --- Obstacle 1 : SPHERE ---
% % % o1.type   = 'sphere';
% % % o1.d_safe = 0.01;
% % % o1.point  = [0.18   0.08   0.06];
% % % o1.radius = 0.04;
% % % obs{end+1} = o1;
% % % 
% % % % --- Obstacle 2 : SPHERE ---
% % % o2.type   = 'sphere';
% % % o2.d_safe = 0.01;
% % % o2.point  = [0.33   0.04   0.06];
% % % o2.radius = 0.04;
% % % obs{end+1} = o2;
% % % 
% % % % --- Obstacle 3 : PLANE ---
% % % o3.type   = 'plane';
% % % o3.d_safe = 0.02;
% % % o3.normal = [0  1  1];
% % % o3.point  = [0 -0.05  0];
% % % obs{end+1} = o3;
% % % 
% % % % % --- Optional: CYLINDER ---
% % % % o4.type   = 'cylinder';
% % % % o4.d_safe = 0.02;
% % % % o4.point  = [0.22   0.05   0.0];
% % % % o4.axis   = [0  0  1];
% % % % o4.radius = 0.03;
% % % % obs{end+1} = o4;
% % % 
% % % % % --- Optional: ELLIPSOID ---
% % % % o5.type     = 'ellipsoid';
% % % % o5.d_safe   = 0.02;
% % % % o5.point    = [0.25   0.12   0.10];
% % % % o5.semiaxes = [0.04   0.03   0.03];
% % % % o5.axes     = eye(3);
% % % % obs{end+1}  = o5;
% % % 
% % % %% -----------------------------------------------------------------------
% % % %  BARRIER + AUGMENTED LAGRANGIAN PARAMETERS
% % % % -----------------------------------------------------------------------
% % % mu_obs   = 5e-3;   % log-barrier weight (shared across obstacles)
% % % mu_tip   = 10;     % initial AL penalty weight
% % % mu_scale = 5;      % AL penalty growth factor per outer iteration
% % % lambda   = zeros(3,1);
% % % nAL      = 12;     % max AL outer iterations
% % % tolTip   = 1e-5;   % convergence tolerance on tip error
% % % 
% % % %% -----------------------------------------------------------------------
% % % %  DISCRETISATION
% % % % -----------------------------------------------------------------------
% % % nPts = 80 * ones(1,m);
% % % 
% % % %% -----------------------------------------------------------------------
% % % %  INITIAL GUESS
% % % % -----------------------------------------------------------------------
% % % x0 = [1.5*ones(1,m)   0.3*ones(1,m)];
% % % x  = x0;
% % % 
% % % %% -----------------------------------------------------------------------
% % % %  OPTIMISER OPTIONS
% % % % -----------------------------------------------------------------------
% % % opts = optimoptions('fminunc', ...
% % %     'Algorithm',              'quasi-newton', ...
% % %     'Display',                'iter',         ...
% % %     'MaxIterations',          600,            ...
% % %     'MaxFunctionEvaluations', 3e5,            ...
% % %     'OptimalityTolerance',    1e-8,           ...
% % %     'StepTolerance',          1e-10);
% % % 
% % % %% -----------------------------------------------------------------------
% % % %  AUGMENTED LAGRANGIAN OUTER LOOP
% % % % -----------------------------------------------------------------------
% % % fprintf('=====================================================\n');
% % % fprintf(' TDCR + Multi-Obstacle Log Barrier + Augmented Lagrangian\n');
% % % fprintf('=====================================================\n\n');
% % % 
% % % for al = 1:nAL
% % % 
% % %     fprintf('AL ITERATION %d\n', al);
% % % 
% % %     obj = @(xx) totalCost(xx, Ls, nPts, EI, GJ, obs, r_des, mu_obs, lambda, mu_tip);
% % % 
% % %     x = fminunc(obj, x, opts);
% % % 
% % %     [~, rline] = rodFK(x(1:m), x(m+1:end), Ls, nPts);
% % %     e = rline(end,:)' - r_des(:);
% % % 
% % %     fprintf('Tip Error = %.6e\n', norm(e));
% % % 
% % %     lambda = lambda + mu_tip * e;
% % %     mu_tip = mu_tip * mu_scale;
% % % 
% % %     fprintf('mu_tip    = %.3e\n', mu_tip);
% % %     fprintf('--------------------------------------------\n');
% % % 
% % %     if norm(e) < tolTip
% % %         fprintf('Converged.\n');
% % %         break;
% % %     end
% % % end
% % % 
% % % %% -----------------------------------------------------------------------
% % % %  FINAL RESULT
% % % % -----------------------------------------------------------------------
% % % kappa = x(1:m);
% % % phi   = x(m+1:end);
% % % 
% % % [~, rline, r_junc] = rodFK(kappa, phi, Ls, nPts);
% % % 
% % % Eel = elasticEnergy(kappa, phi, Ls, EI, GJ);
% % % 
% % % fprintf('\n=====================================================\n');
% % % fprintf('FINAL RESULT\n');
% % % fprintf('=====================================================\n');
% % % fprintf('Elastic Energy    = %.6f\n', Eel);
% % % fprintf('Tip Error         = %.6e\n', norm(rline(end,:) - r_des));
% % % 
% % % fprintf('\nObstacle clearances (positive = safe):\n');
% % % for p = 1:numel(obs)
% % %     d_min_p = minDistToSurfaceHard(rline, obs{p});
% % %     fprintf('  Obs %d (%s):  d_min=%.4f m   d_safe=%.4f m   margin=%.4f m\n', ...
% % %             p, obs{p}.type, d_min_p, obs{p}.d_safe, d_min_p - obs{p}.d_safe);
% % % end
% % % 
% % % %% -----------------------------------------------------------------------
% % % %  VISUALISATION
% % % % -----------------------------------------------------------------------
% % % nPlot = 200 * ones(1,m);
% % % [~, rLine, rJunc] = rodFK(kappa, phi, Ls, nPlot);
% % % 
% % % figure('Color','w','Position',[100 100 900 700]);
% % % hold on; view(3);
% % % 
% % % colors = lines(numel(obs));
% % % for p = 1:numel(obs)
% % %     plotObstacle(obs{p}, colors(p,:));
% % % end
% % % 
% % % h1 = plot3(rLine(:,2), rLine(:,3), rLine(:,1), 'r-', 'LineWidth', 3, ...
% % %       'DisplayName','Rod centreline');
% % % h2 = plot3(rJunc(:,2), rJunc(:,3), rJunc(:,1), '.g', 'MarkerSize', 22, ...
% % %       'DisplayName','Section junctions');
% % % h3 = plot3(r_des(2), r_des(3), r_des(1), '*r', 'MarkerSize', 15, ...
% % %       'LineWidth', 2.5, 'DisplayName', ...
% % %       sprintf('Target tip (err=%.2emm)', norm(rLine(end,:)-r_des)*1e3));
% % % 
% % % grid on; axis equal; box on;
% % % xlabel('Y (m)'); ylabel('Z (m)'); zlabel('X (m)');
% % % zlim([0 0.4]);
% % % set(gca,'FontSize',12);
% % % legend([h1 h2 h3],'Location','bestoutside','Interpreter','latex');
% % % title('TDCR: Log-Barrier + Augmented-Lagrangian Obstacle Avoidance', ...
% % %       'Interpreter','latex','FontSize',14);
% % % camlight; lighting gouraud;
% % % 
% % % %% ========================================================================
% % % %  TOTAL COST  (objective passed to fminunc)
% % % %% ========================================================================
% % % function F = totalCost(x, Ls, nPts, EI, GJ, obs, r_des, mu_obs, lambda, mu_tip)
% % % 
% % % m     = numel(Ls);
% % % kappa = x(1:m);
% % % phi   = x(m+1:end);
% % % 
% % % Eel = elasticEnergy(kappa, phi, Ls, EI, GJ);
% % % 
% % % [~, rline] = rodFK(kappa, phi, Ls, nPts);
% % % 
% % % Eobs = 0;
% % % for p = 1:numel(obs)
% % %     Eobs = Eobs + obstacleBarrier(rline, obs{p}, mu_obs);
% % % end
% % % 
% % % e    = rline(end,:)' - r_des(:);
% % % Etip = lambda' * e + 0.5 * mu_tip * (e' * e);
% % % 
% % % F = Eel + Eobs + Etip;
% % % 
% % % fprintf('Eel=%7.4e  Eobs=%7.4e  Etip=%7.4e\r', Eel, Eobs, Etip);
% % % end
% % % 
% % % %% ========================================================================
% % % %  ELASTIC ENERGY
% % % %% ========================================================================
% % % function E = elasticEnergy(kappa, phi, Ls, EI, GJ)
% % % Ebend  = 0.5 * EI * sum(Ls .* kappa.^2);
% % % Etwist = 0.5 * GJ * sum(phi.^2);
% % % E      = Ebend + Etwist;
% % % end
% % % 
% % % %% ========================================================================
% % % %  GENERAL LOG-BARRIER OVER A SINGLE OBSTACLE
% % % %  ---------------------------------------------------------------
% % % %  Feasible region:  signedDist(p,obs) >= d_safe   (outside the safety
% % % %  shell of the obstacle surface, whatever its type).
% % % %
% % % %    gap = signedDist(p,obs) - d_safe
% % % %
% % % %    gap > 0  ->  feasible   ->  -mu * ln( gap / Rref )
% % % %    gap <= 0 ->  violation  ->   1e6 * gap^2   (hard quadratic penalty)
% % % %
% % % %  Rref is a per-obstacle characteristic length used only to
% % % %  non-dimensionalise the log argument (keeps the barrier well scaled
% % % %  regardless of obstacle size).
% % % %% ========================================================================
% % % function Eb = obstacleBarrier(rline, obs, mu)
% % % N    = size(rline,1);
% % % Rref = obstacleRefLength(obs);
% % % Eb   = 0;
% % % for i = 1:N
% % %     p   = rline(i,:);
% % %     d   = signedDist(p, obs);
% % %     gap = d - obs.d_safe;
% % %     if gap <= 0
% % %         Eb = Eb + 1e6 * gap^2;
% % %     else
% % %         Eb = Eb - mu * log(gap / Rref);
% % %     end
% % % end
% % % end
% % % 
% % % function R = obstacleRefLength(obs)
% % % switch lower(obs.type)
% % %     case 'sphere',    R = obs.radius;
% % %     case 'cylinder',  R = obs.radius;
% % %     case 'ellipsoid', R = mean(obs.semiaxes);
% % %     case 'plane',     R = 0.05;   % nominal length scale, plane has no radius
% % %     otherwise,        R = 0.05;
% % % end
% % % end
% % % 
% % % %% ========================================================================
% % % %  MINIMUM (HARD) DISTANCE ALONG THE ROD TO A SURFACE — for reporting only
% % % %% ========================================================================
% % % function d_min = minDistToSurfaceHard(rline, obs)
% % % N = size(rline,1);
% % % d = zeros(N,1);
% % % for k = 1:N
% % %     d(k) = signedDist(rline(k,:), obs);
% % % end
% % % d_min = min(d);
% % % end
% % % 
% % % %% ========================================================================
% % % %  SIGNED DISTANCE FUNCTIONS  (positive = outside / safe side)
% % % %% ========================================================================
% % % function dist = signedDist(p, obs)
% % % switch lower(obs.type)
% % % 
% % %     case 'plane'
% % %         n    = obs.normal(:)' / norm(obs.normal);
% % %         r0   = obs.point(:)';
% % %         dist = dot(n, p - r0);
% % % 
% % %     case 'sphere'
% % %         c    = obs.point(:)';
% % %         dist = norm(p - c) - obs.radius;
% % % 
% % %     case 'cylinder'
% % %         c    = obs.point(:)';
% % %         a    = obs.axis(:)' / norm(obs.axis);
% % %         v    = p - c;
% % %         v_ax = dot(v,a)*a;
% % %         v_r  = v - v_ax;
% % %         dist = norm(v_r) - obs.radius;
% % % 
% % %     case 'ellipsoid'
% % %         c  = obs.point(:)';
% % %         ab = obs.semiaxes(:)';
% % %         Ef = obs.axes;
% % %         pt = p - c;
% % %         u  = [dot(pt,Ef(1,:)), dot(pt,Ef(2,:)), dot(pt,Ef(3,:))];
% % %         F  = sum((u./ab).^2) - 1;
% % %         gF_local = 2*u./(ab.^2);
% % %         gF_world = gF_local(1)*Ef(1,:) + gF_local(2)*Ef(2,:) + gF_local(3)*Ef(3,:);
% % %         gn = norm(gF_world);
% % %         if gn < 1e-12
% % %             dist = 0;
% % %         else
% % %             dist = F / gn;
% % %         end
% % % 
% % %     otherwise
% % %         error('Unknown obstacle type: %s', obs.type);
% % % end
% % % end
% % % 
% % % %% ========================================================================
% % % %  FORWARD KINEMATICS  (Frenet frame, geometric RK4)
% % % %% ========================================================================
% % % function [s_tot, r_tot, r_junc] = rodFK(kappa, phi, Ls, nPts)
% % % 
% % % m = numel(kappa);
% % % 
% % % curv_funs = cellfun(@(k) @(s) k, num2cell(kappa), 'UniformOutput', false);
% % % tors_funs = repmat({@(s) 0}, 1, m);
% % % 
% % % % Rod base points along +X (consistent with the obstacle placement above)
% % % T0 = [1 0 0];
% % % N0 = [0 1 0]*cos(phi(1)) + [0 0 1]*sin(phi(1));
% % % B0 = cross(T0, N0);
% % % 
% % % [s_tot, r_tot, ~, ~, ~] = rodMulti( ...
% % %     curv_funs, tors_funs, Ls, nPts, phi(2:end), ...
% % %     [0 0 0], T0, N0, B0);
% % % 
% % % cumL   = [0 cumsum(Ls)];
% % % r_junc = zeros(numel(cumL),3);
% % % for j = 1:numel(cumL)
% % %     [~, idx]    = min(abs(s_tot - cumL(j)));
% % %     r_junc(j,:) = r_tot(idx,:);
% % % end
% % % end
% % % 
% % % function [s_tot, r_tot, T_tot, N_tot, B_tot] = rodMulti( ...
% % %     curv_funs, tors_funs, Ls, nPts, phis, r0, T0, N0, B0)
% % % 
% % % m = numel(Ls);
% % % s_tot = []; r_tot = []; T_tot = []; N_tot = []; B_tot = [];
% % % r_prev = r0; T_prev = T0; N_prev = N0; B_prev = B0;
% % % cumL = 0;
% % % 
% % % for i = 1:m
% % %     Li = Ls(i);
% % %     ni = nPts(i);
% % % 
% % %     [s_i, r_i, T_i, N_i, B_i] = rodSegment( ...
% % %         curv_funs{i}, tors_funs{i}, Li, ni, r_prev, T_prev, N_prev, B_prev);
% % % 
% % %     s_i = s_i + cumL;
% % % 
% % %     if i == 1
% % %         s_tot = s_i; r_tot = r_i; T_tot = T_i; N_tot = N_i; B_tot = B_i;
% % %     else
% % %         s_tot = [s_tot; s_i(2:end)];        %#ok
% % %         r_tot = [r_tot; r_i(2:end,:)];      %#ok
% % %         T_tot = [T_tot; T_i(2:end,:)];      %#ok
% % %         N_tot = [N_tot; N_i(2:end,:)];      %#ok
% % %         B_tot = [B_tot; B_i(2:end,:)];      %#ok
% % %     end
% % % 
% % %     cumL = cumL + Li;
% % % 
% % %     if i < m
% % %         ph     = phis(i);
% % %         T_prev = T_i(end,:);
% % %         Nend   = N_i(end,:);
% % %         Bend   = B_i(end,:);
% % % 
% % %         N_prev =  cos(ph)*Nend + sin(ph)*Bend;
% % %         B_prev = -sin(ph)*Nend + cos(ph)*Bend;
% % %         r_prev =  r_i(end,:);
% % %     end
% % % end
% % % end
% % % 
% % % function [s, r, T, N, B] = rodSegment(curv_fun, tors_fun, Lseg, nPts, r0, T0, N0, B0)
% % % 
% % % s = linspace(0, Lseg, nPts)';
% % % 
% % % r = zeros(nPts,3); T = zeros(nPts,3);
% % % N = zeros(nPts,3); B = zeros(nPts,3);
% % % 
% % % r(1,:) = r0; T(1,:) = T0; N(1,:) = N0; B(1,:) = B0;
% % % 
% % % for k = 1:nPts-1
% % %     ds  = s(k+1) - s(k);
% % %     kap = curv_fun(s(k));
% % %     tau = tors_fun(s(k));
% % % 
% % %     y0 = [T(k,:)'; N(k,:)'; B(k,:)'];
% % % 
% % %     f = @(y) [ kap*y(4:6); ...
% % %               -kap*y(1:3) + tau*y(7:9); ...
% % %               -tau*y(4:6)];
% % % 
% % %     k1 = f(y0);
% % %     k2 = f(y0 + 0.5*ds*k1);
% % %     k3 = f(y0 + 0.5*ds*k2);
% % %     k4 = f(y0 +     ds*k3);
% % % 
% % %     y1 = y0 + (ds/6)*(k1 + 2*k2 + 2*k3 + k4);
% % % 
% % %     T(k+1,:) = y1(1:3)' / norm(y1(1:3));
% % %     N(k+1,:) = y1(4:6)' / norm(y1(4:6));
% % %     B(k+1,:) = cross(T(k+1,:), N(k+1,:));
% % %     B(k+1,:) = B(k+1,:) / norm(B(k+1,:));
% % %     N(k+1,:) = cross(B(k+1,:), T(k+1,:));
% % % 
% % %     r(k+1,:) = r(k,:) + ds*T(k,:);
% % % end
% % % end
% % % 
% % % %% ========================================================================
% % % %  OBSTACLE VISUALISATION
% % % %% ========================================================================
% % % function plotObstacle(obs, col)
% % % alpha_val = 0.30;
% % % res = 40;
% % % switch lower(obs.type)
% % % 
% % %     case 'plane'
% % %         n=obs.normal(:)'/norm(obs.normal); r0=obs.point(:)';
% % %         tmp=[0 0 1]; if abs(dot(n,tmp))>0.9, tmp=[1 0 0]; end
% % %         e1=cross(n,tmp); e1=e1/norm(e1); e2=cross(n,e1);
% % %         L=0.30;
% % %         [ug,vg]=meshgrid(linspace(-L,L,res));
% % %         Xp=r0(1)+ug*e1(1)+vg*e2(1);
% % %         Yp=r0(2)+ug*e1(2)+vg*e2(2);
% % %         Zp=r0(3)+ug*e1(3)+vg*e2(3);
% % %         surf(Yp,Zp,Xp,'FaceColor',col,'FaceAlpha',alpha_val, ...
% % %              'EdgeColor','none','DisplayName','Plane');
% % % 
% % %     case 'sphere'
% % %         c=obs.point(:)'; R=obs.radius;
% % %         [xs,ys,zs]=sphere(res);
% % %         surf(c(2)+R*ys,c(3)+R*zs,c(1)+R*xs,'FaceColor',col, ...
% % %              'FaceAlpha',alpha_val,'EdgeColor','none','DisplayName','Sphere');
% % % 
% % %     case 'cylinder'
% % %         c=obs.point(:)'; a=obs.axis(:)'/norm(obs.axis); R=obs.radius; H=0.35;
% % %         tmp=[0 0 1]; if abs(dot(a,tmp))>0.9, tmp=[1 0 0]; end
% % %         u1=cross(a,tmp); u1=u1/norm(u1); u2=cross(a,u1);
% % %         th=linspace(0,2*pi,res); hh=linspace(-H,H,res);
% % %         [TH,HH]=meshgrid(th,hh);
% % %         Xc=c(1)+R*cos(TH)*u1(1)+R*sin(TH)*u2(1)+HH*a(1);
% % %         Yc=c(2)+R*cos(TH)*u1(2)+R*sin(TH)*u2(2)+HH*a(2);
% % %         Zc=c(3)+R*cos(TH)*u1(3)+R*sin(TH)*u2(3)+HH*a(3);
% % %         surf(Yc,Zc,Xc,'FaceColor',col,'FaceAlpha',alpha_val, ...
% % %              'EdgeColor','none','DisplayName','Cylinder');
% % % 
% % %     case 'ellipsoid'
% % %         c=obs.point(:)'; ab=obs.semiaxes(:)'; Ef=obs.axes;
% % %         [xs,ys,zs]=ellipsoid(0,0,0,ab(1),ab(2),ab(3),res);
% % %         sz=size(xs); pts=[xs(:),ys(:),zs(:)]*Ef;
% % %         Xel=reshape(pts(:,1),sz)+c(1);
% % %         Yel=reshape(pts(:,2),sz)+c(2);
% % %         Zel=reshape(pts(:,3),sz)+c(3);
% % %         surf(Yel,Zel,Xel,'FaceColor',col,'FaceAlpha',alpha_val, ...
% % %              'EdgeColor','none','DisplayName','Ellipsoid');
% % % end
% % % end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
clc; clear; close all;

% =========================================================================
%  TDCR SHAPE OPTIMISATION — UNIFIED LOG-BARRIER OBSTACLE LIBRARY
%  -----------------------------------------------------------------------
%  One optimisation engine (log barrier + Augmented Lagrangian, solved
%  with fminunc) that can enforce ANY COMBINATION of:
%
%     - SPHERE     obstacle          (avoid: stay outside a ball)
%     - CYLINDER   obstacle / tube   (avoid: stay outside a pole
%                                      OR confine: stay inside a tube —
%                                      this is the "rod threading through
%                                      a cylindrical channel" case)
%     - PLANE      half-space        (avoid: stay on the safe side)
%     - ELLIPSOID  obstacle          (avoid: stay outside)
%
%  Turn any constraint on/off with obs{k}.active = true/false — use one,
%  two, or all of them together. Every obstacle also carries its own
%  .mode ('avoid' or 'confine') so the SAME cylinder definition can model
%  either "don't hit this pole" or "stay inside this pipe".
%
%  Unknowns:  x = [kappa_1 ... kappa_m   phi_1 ... phi_m]
%
%     F(x) = E_elastic + sum_active_p E_barrier(p) + E_tip (Augmented Lagrangian)
%
% =========================================================================

%% -----------------------------------------------------------------------
%  ROD PARAMETERS
% -----------------------------------------------------------------------
Ls    = [0.10 0.10 0.10 0.10 0.10];   % section lengths [m]
m     = numel(Ls);
L_tot = sum(Ls);

fprintf('Total rod length = %.3f m\n\n', L_tot);

%% -----------------------------------------------------------------------
%  MATERIAL (simple scalar stiffnesses, as in the barrier formulation)
% -----------------------------------------------------------------------
EI = 1.0;
GJ = 0.1;

%% -----------------------------------------------------------------------
%  TARGET TIP
% -----------------------------------------------------------------------
r_des = [0.40  0.1  0.1];

fprintf('||r_des|| = %.3f m   Reachable = %s\n\n', ...
        norm(r_des), mat2str(norm(r_des) <= L_tot));

%% -----------------------------------------------------------------------
%  OBSTACLE / CONSTRAINT LIBRARY
%  -----------------------------------------------------------------------
%  Every obstacle struct supports:
%     .type     'sphere' | 'cylinder' | 'plane' | 'ellipsoid'
%     .mode     'avoid'   -> stay OUTSIDE the surface (default)
%               'confine' -> stay INSIDE the surface (workspace/tube)
%     .d_safe   required clearance [m]
%     .active   true/false  -> toggle this constraint on/off
%     .mu       (optional) per-obstacle barrier weight; if omitted the
%               global mu_obs below is used
%
%  Type-specific fields:
%     sphere    : .point (centre), .radius
%     cylinder  : .axis_point, .axis_dir (unit), .radius, and EITHER of:
%                 .height       = H   (finite cylinder from t=0 to t=H,
%                                       measured from axis_point along
%                                       axis_dir — the common case: a
%                                       pole or tube of a given height)
%                 .axial_range  = [tmin tmax]  (general finite window,
%                                       lets the cylinder start/end
%                                       anywhere along the axis, not
%                                       just at t=0)
%                 Omit BOTH for an infinite cylinder (rarely what you want).
%     plane     : .normal, .point
%     ellipsoid : .point, .semiaxes [a b c], .axes (3x3 rows = principal dirs)
% -----------------------------------------------------------------------
obs = {};

% --- Obstacle 1 : CYLINDER used as a CONFINING TUBE (rod threads through it)
oCyl.type        = 'cylinder';
oCyl.mode        = 'confine';
oCyl.d_safe      = 0.003;
oCyl.axis_point  = [0 0 0];
oCyl.axis_dir    = [1 0 0];          % tube runs along +X
oCyl.radius      = 0.05;
oCyl.axial_range = [0  0.2];        % only confine within this stretch
oCyl.active      = true;
oCyl.mu          = 1e-1;             % stiffer weight for the workspace wall
obs{end+1} = oCyl;

% --- Obstacle 2 : SPHERE obstacle to AVOID (sits inside the tube)
oSph.type   = 'sphere';
oSph.mode   = 'avoid';
oSph.d_safe = 0.005;
oSph.point  = [0.25  0.01  0.01];
oSph.radius = 0.015;
oSph.active = true;
obs{end+1} = oSph;

% --- Obstacle 3 : PLANE to AVOID (e.g. a floor / bulkhead)
oPln.type   = 'plane';
oPln.mode   = 'avoid';
oPln.d_safe = 0.005;
oPln.normal = [0  0  1];             % feasible side is +Z
oPln.point  = [0  0 -0.04];
oPln.active = true;
obs{end+1} = oPln;

% --- Obstacle 4 : CYLINDER pole to AVOID, with a FINITE HEIGHT
%     (a real post/column of height 0.12 m standing up from axis_point,
%     not an infinite rod — this is what .height controls)
% % oPole.type       = 'cylinder';
% % oPole.mode       = 'avoid';
% % oPole.d_safe     = 0.01;
% % oPole.axis_point = [0.30  0.10  0];
% % oPole.axis_dir   = [0  0  1];
% % oPole.radius     = 0.02;
% % oPole.height     = 0.12;      % pole exists only for t in [0, 0.12] along axis_dir
% % oPole.active     = true;
% % obs{end+1} = oPole;

% % --- Optional : an ELLIPSOID obstacle to AVOID
% oEll.type     = 'ellipsoid';
% oEll.mode     = 'avoid';
% oEll.d_safe   = 0.01;
% oEll.point    = [0.3  0.03  -0.01];
% oEll.semiaxes = [0.03  0.02  0.02];
% oEll.axes     = eye(3);
% oEll.active   = false;
% obs{end+1} = oEll;

%% -----------------------------------------------------------------------
%  GLOBAL BARRIER + AUGMENTED LAGRANGIAN PARAMETERS
% -----------------------------------------------------------------------
mu_obs   = 5e-3;   % default log-barrier weight (used unless obs{k}.mu given)
mu_tip   = 10;     % initial AL penalty weight
mu_scale = 5;      % AL penalty growth factor per outer iteration
lambda   = zeros(3,1);
nAL      = 10;
tolTip   = 1e-5;

%% -----------------------------------------------------------------------
%  DISCRETISATION
% -----------------------------------------------------------------------
nPts = 80 * ones(1,m);

%% -----------------------------------------------------------------------
%  INITIAL GUESS
% -----------------------------------------------------------------------
x0 = [1.5*ones(1,m)   zeros(1,m)];
x  = x0;

%% -----------------------------------------------------------------------
%  OPTIMISER OPTIONS
% -----------------------------------------------------------------------
opts = optimoptions('fminunc', ...
    'Algorithm',              'quasi-newton', ...
    'Display',                'iter',         ...
    'MaxIterations',          600,            ...
    'MaxFunctionEvaluations', 3e5,            ...
    'OptimalityTolerance',    1e-8,           ...
    'StepTolerance',          1e-10);

%% -----------------------------------------------------------------------
%  ACTIVE-OBSTACLE FILTER  — this is the "choose one, two, or all" switch
% -----------------------------------------------------------------------
activeMask = cellfun(@(o) isfield(o,'active') && o.active, obs);
obsActive  = obs(activeMask);

fprintf('Active constraints (%d of %d):\n', numel(obsActive), numel(obs));
for p = 1:numel(obsActive)
    fprintf('  - %s (%s)\n', obsActive{p}.type, obsActive{p}.mode);
end
fprintf('\n');

%% -----------------------------------------------------------------------
%  AUGMENTED LAGRANGIAN OUTER LOOP
% -----------------------------------------------------------------------
fprintf('=====================================================\n');
fprintf(' TDCR + Unified Log-Barrier Library + Aug. Lagrangian\n');
fprintf('=====================================================\n\n');

for al = 1:nAL

    fprintf('AL ITERATION %d\n', al);

    obj = @(xx) totalCost(xx, Ls, nPts, EI, GJ, obsActive, r_des, mu_obs, lambda, mu_tip);

    x = fminunc(obj, x, opts);

    [~, rline] = rodFK(x(1:m), x(m+1:end), Ls, nPts);
    e = rline(end,:)' - r_des(:);

    fprintf('Tip Error = %.6e\n', norm(e));

    lambda = lambda + mu_tip * e;
    mu_tip = mu_tip * mu_scale;

    fprintf('mu_tip    = %.3e\n', mu_tip);
    fprintf('--------------------------------------------\n');

    if norm(e) < tolTip
        fprintf('Converged.\n');
        break;
    end
end

%% -----------------------------------------------------------------------
%  FINAL RESULT
% -----------------------------------------------------------------------
kappa = x(1:m);
phi   = x(m+1:end);

[~, rline, r_junc] = rodFK(kappa, phi, Ls, nPts);
Eel = elasticEnergy(kappa, phi, Ls, EI, GJ);

fprintf('\n=====================================================\n');
fprintf('FINAL RESULT\n');
fprintf('=====================================================\n');
fprintf('Elastic Energy = %.6f\n', Eel);
fprintf('Tip Error      = %.6e\n', norm(rline(end,:) - r_des));

fprintf('\nConstraint margins (positive = satisfied):\n');
for p = 1:numel(obsActive)
    m_p = worstMarginHard(rline, obsActive{p});
    fprintf('  %s (%s): worst-point margin = %.4f m  (d_safe = %.4f m)\n', ...
            obsActive{p}.type, obsActive{p}.mode, m_p, obsActive{p}.d_safe);
end

%% -----------------------------------------------------------------------
%  VISUALISATION
% -----------------------------------------------------------------------
nPlot = 200 * ones(1,m);
[~, rLine, rJunc] = rodFK(kappa, phi, Ls, nPlot);

figure('Color','w','Position',[100 100 900 700]);
hold on; view(3);

colors = lines(numel(obsActive));
for p = 1:numel(obsActive)
    plotObstacle(obsActive{p}, colors(p,:));
end

h1 = plot3(rLine(:,1), rLine(:,2), rLine(:,3), 'r-', 'LineWidth', 3, ...
      'DisplayName','Rod centreline');
h2 = plot3(rJunc(:,1), rJunc(:,2), rJunc(:,3), '.g', 'MarkerSize', 22, ...
      'DisplayName','Section junctions');
h3 = plot3(r_des(1), r_des(2), r_des(3), '*r', 'MarkerSize', 15, ...
      'LineWidth', 2.5, 'DisplayName', ...
      sprintf('Target tip (err=%.2emm)', norm(rLine(end,:)-r_des)*1e3));

grid on; axis equal; box on;
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
set(gca,'FontSize',12);
legend([h1 h2 h3],'Location','bestoutside','Interpreter','latex');
title('TDCR: Unified Log-Barrier Obstacle Library (Sphere / Cylinder / Plane)', ...
      'Interpreter','latex','FontSize',14);
camlight; lighting gouraud;

%% ========================================================================
%  TOTAL COST
%% ========================================================================
function F = totalCost(x, Ls, nPts, EI, GJ, obs, r_des, mu_obs_default, lambda, mu_tip)

m     = numel(Ls);
kappa = x(1:m);
phi   = x(m+1:end);

Eel = elasticEnergy(kappa, phi, Ls, EI, GJ);

[~, rline] = rodFK(kappa, phi, Ls, nPts);

Eobs = 0;
for p = 1:numel(obs)
    mu_p = mu_obs_default;
    if isfield(obs{p},'mu') && ~isempty(obs{p}.mu)
        mu_p = obs{p}.mu;
    end
    Eobs = Eobs + obstacleBarrier(rline, obs{p}, mu_p);
end

e    = rline(end,:)' - r_des(:);
Etip = lambda' * e + 0.5 * mu_tip * (e' * e);

F = Eel + Eobs + Etip;

fprintf('Eel=%7.4e  Eobs=%7.4e  Etip=%7.4e\r', Eel, Eobs, Etip);
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
%  GENERAL LOG-BARRIER OVER A SINGLE OBSTACLE (avoid OR confine)
%  ---------------------------------------------------------------
%  For every rod point p:
%     [d, active_here] = surfaceDistance(p, obs)
%        d > 0 always means "on the safe side of the raw surface", i.e.
%        outside for spheres/cylinders/ellipsoids/planes as normally
%        defined. active_here is false only for cylinder points that
%        fall outside an optional .axial_range window (no constraint
%        applies there).
%
%     mode = 'avoid'   -> gap =  d - d_safe   (must stay outside)
%     mode = 'confine' -> gap = -d - d_safe   (must stay inside)
%
%     gap > 0  -> feasible -> -mu * log(gap / Rref)
%     gap <= 0 -> violation -> 1e6 * gap^2
%% ========================================================================
function Eb = obstacleBarrier(rline, obs, mu)
N    = size(rline,1);
Rref = obstacleRefLength(obs);
mode = 'avoid';
if isfield(obs,'mode') && ~isempty(obs.mode)
    mode = lower(obs.mode);
end

Eb = 0;
for i = 1:N
    p = rline(i,:);
    [d, isActiveHere] = surfaceDistance(p, obs);
    if ~isActiveHere
        continue;   % point outside the obstacle's axial window: no constraint
    end

    if strcmp(mode,'confine')
        gap = -d - obs.d_safe;
    else
        gap =  d - obs.d_safe;
    end

    if gap <= 0
        Eb = Eb + 1e6 * gap^2;
    else
        Eb = Eb - mu * log(gap / Rref);
    end
end
end

function R = obstacleRefLength(obs)
switch lower(obs.type)
    case 'sphere',    R = obs.radius;
    case 'cylinder',  R = obs.radius;
    case 'ellipsoid', R = mean(obs.semiaxes);
    case 'plane',     R = 0.05;   % nominal scale, plane has no radius
    otherwise,        R = 0.05;
end
end

%% ========================================================================
%  WORST-CASE (HARD) MARGIN ALONG THE ROD — for reporting only
%% ========================================================================
function m_worst = worstMarginHard(rline, obs)
N    = size(rline,1);
mode = 'avoid';
if isfield(obs,'mode') && ~isempty(obs.mode)
    mode = lower(obs.mode);
end
m_worst = inf;
for k = 1:N
    [d, isActiveHere] = surfaceDistance(rline(k,:), obs);
    if ~isActiveHere
        continue;
    end
    if strcmp(mode,'confine')
        gap = -d - obs.d_safe;
    else
        gap =  d - obs.d_safe;
    end
    m_worst = min(m_worst, gap);
end
if isinf(m_worst)
    m_worst = NaN;   % obstacle never applied along this rod (e.g. axial window missed entirely)
end
end

%% ========================================================================
%  RAW SURFACE DISTANCE  (positive = outside, in the natural sense of the
%  surface itself — mode handling happens in obstacleBarrier)
%  Second output flags whether the constraint applies at this point
%  (only relevant for cylinders with a finite .axial_range).
%% ========================================================================
function [dist, isActiveHere] = surfaceDistance(p, obs)
isActiveHere = true;

switch lower(obs.type)

    case 'plane'
        n    = obs.normal(:)' / norm(obs.normal);
        r0   = obs.point(:)';
        dist = dot(n, p - r0);

    case 'sphere'
        c    = obs.point(:)';
        dist = norm(p - c) - obs.radius;

    case 'cylinder'
        c = obs.axis_point(:)';
        a = obs.axis_dir(:)' / norm(obs.axis_dir);
        v = p - c;
        t = dot(v,a);

        [tmin, tmax, isBounded] = cylinderAxialBounds(obs);
        if isBounded && (t < tmin || t > tmax)
            dist = 0;
            isActiveHere = false;
            return;
        end

        v_r  = v - t*a;
        dist = norm(v_r) - obs.radius;

    case 'ellipsoid'
        c  = obs.point(:)';
        ab = obs.semiaxes(:)';
        Ef = obs.axes;
        pt = p - c;
        u  = [dot(pt,Ef(1,:)), dot(pt,Ef(2,:)), dot(pt,Ef(3,:))];
        F  = sum((u./ab).^2) - 1;
        gF_local = 2*u./(ab.^2);
        gF_world = gF_local(1)*Ef(1,:) + gF_local(2)*Ef(2,:) + gF_local(3)*Ef(3,:);
        gn = norm(gF_world);
        if gn < 1e-12
            dist = 0;
        else
            dist = F / gn;
        end

    otherwise
        error('Unknown obstacle type: %s', obs.type);
end
end

%% ========================================================================
%  RESOLVE A CYLINDER'S AXIAL EXTENT
%  Priority: explicit .axial_range  >  .height (treated as [0, height])  >
%  unbounded (infinite cylinder).
%% ========================================================================
function [tmin, tmax, isBounded] = cylinderAxialBounds(obs)
if isfield(obs,'axial_range') && ~isempty(obs.axial_range)
    tmin = obs.axial_range(1);
    tmax = obs.axial_range(2);
    isBounded = true;
elseif isfield(obs,'height') && ~isempty(obs.height)
    tmin = 0;
    tmax = obs.height;
    isBounded = true;
else
    tmin = -inf; tmax = inf;
    isBounded = false;
end
end

%% ========================================================================
%  FORWARD KINEMATICS
%% ========================================================================
function [s_tot, r_tot, r_junc] = rodFK(kappa, phi, Ls, nPts)

m = numel(kappa);

curv_funs = cellfun(@(k) @(s) k, num2cell(kappa), 'UniformOutput', false);
tors_funs = repmat({@(s) 0}, 1, m);

T0 = [1 0 0];                              % rod grows along +X
N0 = [0 1 0]*cos(phi(1)) + [0 0 1]*sin(phi(1));
B0 = cross(T0, N0);

[s_tot, r_tot, ~, ~, ~] = rodMulti( ...
    curv_funs, tors_funs, Ls, nPts, phi(2:end), ...
    [0 0 0], T0, N0, B0);

cumL   = [0 cumsum(Ls)];
r_junc = zeros(numel(cumL),3);
for j = 1:numel(cumL)
    [~, idx]    = min(abs(s_tot - cumL(j)));
    r_junc(j,:) = r_tot(idx,:);
end
end

function [s_tot, r_tot, T_tot, N_tot, B_tot] = rodMulti( ...
    curv_funs, tors_funs, Ls, nPts, phis, r0, T0, N0, B0)

m = numel(Ls);
s_tot = []; r_tot = []; T_tot = []; N_tot = []; B_tot = [];
r_prev = r0; T_prev = T0; N_prev = N0; B_prev = B0;
cumL = 0;

for i = 1:m
    Li = Ls(i);
    ni = nPts(i);

    [s_i, r_i, T_i, N_i, B_i] = rodSegment( ...
        curv_funs{i}, tors_funs{i}, Li, ni, r_prev, T_prev, N_prev, B_prev);

    s_i = s_i + cumL;

    if i == 1
        s_tot = s_i; r_tot = r_i; T_tot = T_i; N_tot = N_i; B_tot = B_i;
    else
        s_tot = [s_tot; s_i(2:end)];        %#ok
        r_tot = [r_tot; r_i(2:end,:)];      %#ok
        T_tot = [T_tot; T_i(2:end,:)];      %#ok
        N_tot = [N_tot; N_i(2:end,:)];      %#ok
        B_tot = [B_tot; B_i(2:end,:)];      %#ok
    end

    cumL = cumL + Li;

    if i < m
        ph     = phis(i);
        T_prev = T_i(end,:);
        Nend   = N_i(end,:);
        Bend   = B_i(end,:);

        N_prev =  cos(ph)*Nend + sin(ph)*Bend;
        B_prev = -sin(ph)*Nend + cos(ph)*Bend;
        r_prev =  r_i(end,:);
    end
end
end

function [s, r, T, N, B] = rodSegment(curv_fun, tors_fun, Lseg, nPts, r0, T0, N0, B0)

s = linspace(0, Lseg, nPts)';

r = zeros(nPts,3); T = zeros(nPts,3);
N = zeros(nPts,3); B = zeros(nPts,3);

r(1,:) = r0; T(1,:) = T0; N(1,:) = N0; B(1,:) = B0;

for k = 1:nPts-1
    ds  = s(k+1) - s(k);
    kap = curv_fun(s(k));
    tau = tors_fun(s(k));

    y0 = [T(k,:)'; N(k,:)'; B(k,:)'];

    f = @(y) [ kap*y(4:6); ...
              -kap*y(1:3) + tau*y(7:9); ...
              -tau*y(4:6)];

    k1 = f(y0);
    k2 = f(y0 + 0.5*ds*k1);
    k3 = f(y0 + 0.5*ds*k2);
    k4 = f(y0 +     ds*k3);

    y1 = y0 + (ds/6)*(k1 + 2*k2 + 2*k3 + k4);

    T(k+1,:) = y1(1:3)' / norm(y1(1:3));
    N(k+1,:) = y1(4:6)' / norm(y1(4:6));
    B(k+1,:) = cross(T(k+1,:), N(k+1,:));
    B(k+1,:) = B(k+1,:) / norm(B(k+1,:));
    N(k+1,:) = cross(B(k+1,:), T(k+1,:));

    r(k+1,:) = r(k,:) + ds*T(k,:);
end
end

%% ========================================================================
%  OBSTACLE VISUALISATION  (handles 'avoid' and 'confine' cylinders,
%  finite axial_range if given, and translucent solid shells for all types)
%% ========================================================================
function plotObstacle(obs, col)
alpha_val = 0.25;
res = 40;

switch lower(obs.type)

    case 'plane'
        n=obs.normal(:)'/norm(obs.normal); r0=obs.point(:)';
        tmp=[0 0 1]; if abs(dot(n,tmp))>0.9, tmp=[1 0 0]; end
        e1=cross(n,tmp); e1=e1/norm(e1); e2=cross(n,e1);
        L=0.30;
        [ug,vg]=meshgrid(linspace(-L,L,res));
        Xp=r0(1)+ug*e1(1)+vg*e2(1);
        Yp=r0(2)+ug*e1(2)+vg*e2(2);
        Zp=r0(3)+ug*e1(3)+vg*e2(3);
        surf(Xp,Yp,Zp,'FaceColor',col,'FaceAlpha',alpha_val, ...
             'EdgeColor','none','DisplayName','Plane');

    case 'sphere'
        c=obs.point(:)'; R=obs.radius;
        [xs,ys,zs]=sphere(res);
        surf(c(1)+R*xs,c(2)+R*ys,c(3)+R*zs,'FaceColor',col, ...
             'FaceAlpha',alpha_val,'EdgeColor','none','DisplayName','Sphere');

    case 'cylinder'
        c = obs.axis_point(:)';
        a = obs.axis_dir(:)'/norm(obs.axis_dir);
        R = obs.radius;

        [tmin, tmax, isBounded] = cylinderAxialBounds(obs);
        if ~isBounded
            tmin = -0.35; tmax = 0.35;   % just a display extent for a truly infinite cylinder
        end

        tmp=[0 0 1]; if abs(dot(a,tmp))>0.9, tmp=[1 0 0]; end
        u1=cross(a,tmp); u1=u1/norm(u1); u2=cross(a,u1);
        th=linspace(0,2*pi,res); hh=linspace(tmin,tmax,res);
        [TH,HH]=meshgrid(th,hh);
        Xc=c(1)+R*cos(TH)*u1(1)+R*sin(TH)*u2(1)+HH*a(1);
        Yc=c(2)+R*cos(TH)*u1(2)+R*sin(TH)*u2(2)+HH*a(2);
        Zc=c(3)+R*cos(TH)*u1(3)+R*sin(TH)*u2(3)+HH*a(3);

        isConfine = isfield(obs,'mode') && strcmpi(obs.mode,'confine');
        if isConfine
            % draw as a translucent tube (workspace boundary)
            surf(Xc,Yc,Zc,'FaceColor',col,'FaceAlpha',0.12, ...
                 'EdgeColor',col*0.6,'EdgeAlpha',0.25,'LineStyle',':', ...
                 'DisplayName','Cylinder (confine)');
        else
            surf(Xc,Yc,Zc,'FaceColor',col,'FaceAlpha',alpha_val, ...
                 'EdgeColor','none','DisplayName','Cylinder (avoid)');
        end

    case 'ellipsoid'
        c=obs.point(:)'; ab=obs.semiaxes(:)'; Ef=obs.axes;
        [xs,ys,zs]=ellipsoid(0,0,0,ab(1),ab(2),ab(3),res);
        sz=size(xs); pts=[xs(:),ys(:),zs(:)]*Ef;
        Xel=reshape(pts(:,1),sz)+c(1);
        Yel=reshape(pts(:,2),sz)+c(2);
        Zel=reshape(pts(:,3),sz)+c(3);
        surf(Xel,Yel,Zel,'FaceColor',col,'FaceAlpha',alpha_val, ...
             'EdgeColor','none','DisplayName','Ellipsoid');
end
end