%% Problem 1a: Barycenter N Body Problem
% N-Body Solar System Simulation (AU, years, solar masses)

clear; clc; close all;

% Constants
G = 4*pi^2;   % AU^3 / (yr^2 * solar mass)

% Masses (Sun + planets)
mass_bodies = [
    1.0;          % Sun
    3.003e-6;     % Earth
    3.213e-7;     % Mars
    9.545e-4;     % Jupiter
];

% Initial Positions (AU) 
r0 = [
    0; 0; 0;        % Sun
    1; 0; 0;        % Earth
    1.524; 0; 0;    % Mars
    5.204; 0; 0;    % Jupiter
];

% Initial Velocities (AU/year) 
v0 = [
    0; 0; 0;                        % Sun
    0; 2*pi; 0;                     % Earth
    0; 2*pi/sqrt(1.524); 0;         % Mars
    0; 2*pi/sqrt(5.204); 0;         % Jupiter
];

n = length(mass_bodies);
M = sum(mass_bodies);

% Reshape into 3 x n
r0_reshape = reshape(r0,3,n);
v0_reshape = reshape(v0,3,n);

% Compute COM
R = (r0_reshape * mass_bodies(:)) / M;
V = (v0_reshape * mass_bodies(:)) / M;

% Shift to barycentric frame
script_r0 = r0_reshape - R;
script_v0 = v0_reshape - V;

% Now build state vector
X0 = [script_r0(:); script_v0(:)];

% Time span (years)
tspan = [0 20];

% ODE Solver options for stability
options = odeset('RelTol',1e-10,'AbsTol',1e-12);

% Run N-body problem eq. of motion through solver
[t, X] = ode45(@(t,X) Eq_of_motion_NBody(t,X,mass_bodies,G), tspan, X0, options);

% 3D Trajectory Plot
names = {'Sun','Earth','Mars','Jupiter','N+1st Body, Satellite'};
figure; hold on;

colors = lines(n);

% Legend handles (only for marker meaning)
h_init = plot3(nan, nan, nan, 'o', ...
    'Color', 'w', 'MarkerFaceColor', 'w', 'MarkerSize', 7);

h_final = plot3(nan, nan, nan, 's', ...
    'Color', 'w', 'MarkerFaceColor', 'w', 'MarkerSize', 7);

for i = 1:n
    
    x_i = X(:,3*(i-1)+1);
    y_i = X(:,3*(i-1)+2);
    z_i = X(:,3*(i-1)+3);
    
    % Trajectory
    plot3(x_i, y_i, z_i, 'Color', colors(i,:), 'LineWidth', 1.5);
    
    % Initial position
    plot3(x_i(1), y_i(1), z_i(1), 'o', ...
        'Color', colors(i,:), 'MarkerFaceColor', colors(i,:));
    
    % Final position
    plot3(x_i(end), y_i(end), z_i(end), 's', ...
        'Color', colors(i,:), 'MarkerFaceColor', colors(i,:));
    
    % Body label (clean + readable)
    text(x_i(end), y_i(end), z_i(end), ...
        ['  ' names{i}], ...
        'Color', colors(i,:), ...
        'FontSize', 10, ...
        'FontWeight', 'bold');
end

grid on
axis equal
xlabel('x (AU)')
ylabel('y (AU)')
zlabel('z (AU)')
title('3D N-Body Barycentric Trajectories (AU, Years, Solar Masses)')

legend([h_init, h_final], ...
    {'Initial Position','Final Position'})

view(3)

% Check angular momentum

[H,max_rel_drift] = check_angular_momentum(X,t,mass_bodies);

% Plot Angular Momentum Components over time
figure; hold on;

plot(t, H(:,1), 'r', 'LineWidth', 1.5)   % x component
plot(t, H(:,2), 'g--', 'LineWidth', 1.5) % y component 
plot(t, H(:,3), 'b', 'LineWidth', 1.5)   % z component

xlabel('Time (years)')
ylabel('Angular Momentum (AU^2 * solar mass / year)')
title('Angular Momentum Components vs Time')

legend('H_x','H_y','H_z')
grid on

% Check numerical error
disp(['Max relative angular momentum drift: ', num2str(max_rel_drift)])

% Check Energy
[T, U, E] = check_energy(t, X, n, mass_bodies, G); 

% Relative Energy Error plot

% Calculate relative error in energy over time
rel_error = abs(E - E(1)) / abs(E(1));

figure;
semilogy(t, rel_error, 'LineWidth', 2)

xlabel('Time (years)')
ylabel('|E(t) - E(0)| / |E(0)|')
title('Relative Energy Error (Log Scale)')
grid on

% The barycentric N-body integrator conserves both total energy and angular 
% momentum to machine precision, with a maximum relative angular momentum 
% drift of 8.6e-13. This demonstrates that the numerical integration is 
% accurate and that the implementation of the gravitational interactions is 
% correct. The trajectories and invariants confirm the expected physical 
% behavior of an isolated N-body system in an inertial (COM) frame.

%% 1.b N+1st restricted NBP

% Pre-computed NBP
t_NBP_bary = t;
X_NBP_bary = X;

n = length(mass_bodies);

% Build interpolants for barycentric body trajectories History
r_interp = cell(n,1);

for j = 1:n
    
    % Position history of body j
    rj_hist = X_NBP_bary(:,3*(j-1)+1:3*j);

    % Create interpolant
    r_interp{j} = griddedInterpolant(t_NBP_bary, rj_hist, ...
        'linear', 'nearest');

end

% Initial spacecraft state in barycentric frame
scriptR0_sat_bary = [5;0;0] - R;
scriptV0_sat_bary = [0;2*pi/sqrt(50);0] - V;
X0_nplus1_bary = [scriptR0_sat_bary(:); scriptV0_sat_bary(:)];

% Run N+1st-barycentric problem eq. of motion through solver
[tq_bary, Xp_bary] = ode45(@(t,X) ...
    restricted_NBP_barycentric(t, X, mass_bodies, r_interp, G), ...
    tspan, X0_nplus1_bary, options);

% 3D Trajectory Plot
names = {'Sun','Earth','Mars','Jupiter','N+1st Body: Satellite'};
figure; hold on;

colors = lines(n);

% Legend handles (only for marker meaning)
h_init = plot3(nan, nan, nan, 'o', ...
    'Color', 'w', 'MarkerFaceColor', 'w', 'MarkerSize', 7);

h_final = plot3(nan, nan, nan, 's', ...
    'Color', 'w', 'MarkerFaceColor', 'w', 'MarkerSize', 7);

% Plot particle
x_p = Xp_bary(:,1);
y_p = Xp_bary(:,2);
z_p = Xp_bary(:,3);
h_nplus1 = plot3(x_p, y_p, z_p, 'w', 'LineWidth', 2);
plot3(x_p(1), y_p(1), z_p(1), 'wo', 'MarkerFaceColor', 'w');
plot3(x_p(end), y_p(end), z_p(end), 'ws', 'MarkerFaceColor', 'w');

for i = 1:n
    
    x_i = X(:,3*(i-1)+1);
    y_i = X(:,3*(i-1)+2);
    z_i = X(:,3*(i-1)+3);
    
    % Trajectory
    plot3(x_i, y_i, z_i, 'Color', colors(i,:), 'LineWidth', 1.5);
    
    % Initial position
    plot3(x_i(1), y_i(1), z_i(1), 'o', ...
        'Color', colors(i,:), 'MarkerFaceColor', colors(i,:));
    
    % Final position
    plot3(x_i(end), y_i(end), z_i(end), 's', ...
        'Color', colors(i,:), 'MarkerFaceColor', colors(i,:));
    
    % Body label (clean + readable)
    text(x_i(end), y_i(end), z_i(end), ...
        ['  ' names{i}], ...
        'Color', colors(i,:), ...
        'FontSize', 10, ...
        'FontWeight', 'bold');
end

grid on
axis equal
xlabel('x (AU)')
ylabel('y (AU)')
zlabel('z (AU)')
title('3D N+1st-Body Barycentric Trajectories (AU, Years, Solar Masses)')

legend([h_nplus1, h_init, h_final],{'N+1st Body -- Satellite','Initial Position','Final Position'},'location','southeast')

view(3)

% Check angular momentum

N = size(Xp_bary,1);

H = zeros(N,3);

for k = 1:N
    r = Xp_bary(k,1:3);
    v = Xp_bary(k,4:6);

    H(k,:) = cross(r, v);
end

H0 = H(1,:);
rel_drift = vecnorm(H - H0,2,2) / norm(H0);

max_rel_drift = max(rel_drift);

% Plot Angular Momentum Components over time
figure; hold on;

plot(tq_bary, H(:,1), 'r', 'LineWidth', 1.5)   % x component
plot(tq_bary, H(:,2), 'g--', 'LineWidth', 1.5) % y component 
plot(tq_bary, H(:,3), 'b', 'LineWidth', 1.5)   % z component

xlabel('Time (years)')
ylabel('Angular Momentum (AU^2 * solar mass / year)')
title('Angular Momentum Components vs Time')

legend('H_x','H_y','H_z')
grid on

% Check numerical error
disp(['Max relative angular momentum drift in (N+1)st body: ', num2str(max_rel_drift)])

% Check Energy
T = zeros(length(tq_bary),1);
U = zeros(length(tq_bary),1);
m_p = 1e-27;

for k = 1:length(tq_bary)
    
    % Calculate Kinetic Energy 
    v = Xp_bary(k,4:6); % velocity 
    T(k) = 0.5 * m_p * dot(v, v); % 0.5*mv^2
    
    % Calculate Potential Energy 
    r = Xp_bary(k,1:3); % position of "ith" body
    U_k = 0;
        
    for j = 1:n
        r_j = X(k,3*(j-1)+1 : 3*j); % position of the other 3 bodies
        r_ij = norm(r - r_j); % difference in position
        
        % U = G*mi*mj/rij
        U_k = U_k + m_p* G * mass_bodies(j) / r_ij; 
    end

    U(k) = U_k;
end
   
E = T - U; % Total energy

% Relative Energy Error plot

% Calculate relative error in energy over time
rel_error = abs(E - E(1)) / abs(E(1));

figure;
semilogy(tq_bary, rel_error, 'LineWidth', 2)

xlabel('Time (years)')
ylabel('|E(t) - E(0)| / |E(0)|')
title('Relative Energy Error (Log Scale)')
grid on

% In the restricted (N+1)-body barycentric simulation, the additional body 
% does not influence the motion of the primary N bodies. As a result, the 
% energy and angular momentum of the particle are not conserved, with a 
% relative angular momentum drift of 1.29e-1. This is expected since the 
% system is no longer dynamically closed, and the particle evolves in a 
% time-varying gravitational field.

%% 2.a General Bodycentric Integrator

% Represent initial conditions in bodycentric form
[~,idxN] = max(mass_bodies); % Largest body 'Nth" body
r_N = r0_reshape(:,idxN);
v_N = v0_reshape(:,idxN);

for i = 1:n
    r0_reshape(:,i) = r0_reshape(:,i) - r_N;
    v0_reshape(:,i) = v0_reshape(:,i) - v_N;
end

X0_body = [r0_reshape(:); v0_reshape(:)];

% Run N-body problem eq. of motion through solver
[t, X] = ode45(@(t,X) Eq_of_motion_NBody_bodycentric(t,X,mass_bodies,G), tspan, X0_body, options);

% 3D Trajectory Plot
names = {'Sun (Nth Body)','Earth','Mars','Jupiter'};
figure; hold on;

colors = lines(n);

% Legend handles (only for marker meaning)
h_init = plot3(nan, nan, nan, 'o', ...
    'Color', 'w', 'MarkerFaceColor', 'w', 'MarkerSize', 7);

h_final = plot3(nan, nan, nan, 's', ...
    'Color', 'w', 'MarkerFaceColor', 'w', 'MarkerSize', 7);

for i = 1:n
    
    x_i = X(:,3*(i-1)+1);
    y_i = X(:,3*(i-1)+2);
    z_i = X(:,3*(i-1)+3);
    
    % Trajectory
    plot3(x_i, y_i, z_i, 'Color', colors(i,:), 'LineWidth', 1.5);
    
    % Initial position
    plot3(x_i(1), y_i(1), z_i(1), 'o', ...
        'Color', colors(i,:), 'MarkerFaceColor', colors(i,:));
    
    % Final position
    plot3(x_i(end), y_i(end), z_i(end), 's', ...
        'Color', colors(i,:), 'MarkerFaceColor', colors(i,:));
    
    % Body label (clean + readable)
    text(x_i(end), y_i(end), z_i(end), ...
        ['  ' names{i}], ...
        'Color', colors(i,:), ...
        'FontSize', 10, ...
        'FontWeight', 'bold');
end

grid on
axis equal
xlabel('x (AU)')
ylabel('y (AU)')
zlabel('z (AU)')
title('3D N-Body Bodycentric Trajectories (AU, Years, Solar Masses)')

legend([h_init, h_final], ...
    {'Initial Position','Final Position'})

view(3)

% Compare barycentric → bodycentric for checking integral accuracy 
t_NBP_body = t;
X_NBP_body = X;
Nt = length(t_NBP_body);

% Build interpolants from barycentric solution
r_interp_all = griddedInterpolant(t_NBP_bary, ...
    X_NBP_bary(:,1:3*n), ...
    'linear', 'nearest');

v_interp_all = griddedInterpolant(t_NBP_bary, ...
    X_NBP_bary(:,3*n+1:end), ...
    'linear', 'nearest');

r_interp_N = griddedInterpolant(t_NBP_bary, ...
    X_NBP_bary(:,3*(idxN-1)+1:3*idxN), ...
    'linear', 'nearest');

v_interp_N = griddedInterpolant(t_NBP_bary, ...
    X_NBP_bary(:,3*n + 3*(idxN-1)+1:3*n + 3*idxN), ...
    'linear', 'nearest');

% Preallocate
err = zeros(Nt,1);
err_v = zeros(Nt,1);

for k = 1:Nt
    
    % Get bodycentric positions and velocities at current time
    r_body = X_NBP_body(k,1:3*n);
    v_body = X_NBP_body(k,3*n+1:end);
    
    % Interpolate central body (from barycentric solution) at same time
    r_N = r_interp_N(t_NBP_body(k))';
    v_N = v_interp_N(t_NBP_body(k))';
    
    % Get barycentric positions and velocities at same time
    r_bary = r_interp_all(t_NBP_body(k));
    v_bary = v_interp_all(t_NBP_body(k));
    
    % Shift from barycentric frame → bodycentric frame
    r_true_body = r_bary - repmat(r_N',1,n);
    v_true_body = v_bary - repmat(v_N',1,n);
    
    % Compute position and velocity error between methods
    err(k) = norm(r_true_body - r_body);
    err_v(k) = norm(v_true_body - v_body);
end

% Report max error
disp(['Max relative error between bodycentric and barycentric position: ', num2str(max(err))])
disp(['Max relative error between bodycentric and barycentric velocity: ', num2str(max(err_v))])

% The barycentric solution was transformed  into the bodycentric frame by adding the interpolated central body trajectory. 
% The maximum position error between the two formulations was 9.45×10−6AU, showing that it closely matches the barycentric 
% simulation and confirming the numerical accuracy of the integrator.

%% 2.b N+1st-Body Bodycentric formulation

% Represent initial conditions in bodycentric form
[m_N,idxN] = max(mass_bodies); % Largest body 'Nth" body
r_N = r0_reshape(:,idxN);
v_N = v0_reshape(:,idxN);

% Pre-computed NBP
t_NBP_body = t;
X_NBP_body = X;

n = length(mass_bodies);

% Build interpolants for body trajectories
r_interp = cell(n,1);

for j = 1:n
    
    % Position history of body j
    rj_hist = X_NBP_body(:,3*(j-1)+1:3*j);

    % Create interpolant
    r_interp{j} = griddedInterpolant(t_NBP_body, rj_hist, ...
        'linear', 'nearest');

end

% Initial spacecraft state in bodycentric form
scriptR0_sat_body = [5;0;0] - r_N;
scriptV0_sat_body = [0;2*pi/sqrt(50);0] - v_N;

X0_nplus1_body = [scriptR0_sat_body(:); scriptV0_sat_body(:)];

% Run N+1st-body problem eq. of motion through solver
[tq_body, Xp_body] = ode45(@(t,X) ...
    restricted_NBP_bodycentric(t, X, mass_bodies, r_interp, G), ...
    tspan, X0_nplus1_body, options);

% 3D Trajectory Plot
names = {'Sun (Nth Body)','Earth','Mars','Jupiter','N+1st Body: Satellite'};
figure; hold on;

colors = lines(n);

% Legend handles (only for marker meaning)
h_init = plot3(nan, nan, nan, 'o', ...
    'Color', 'w', 'MarkerFaceColor', 'w', 'MarkerSize', 7);

h_final = plot3(nan, nan, nan, 's', ...
    'Color', 'w', 'MarkerFaceColor', 'w', 'MarkerSize', 7);

% Plot particle
x_p = Xp_body(:,1);
y_p = Xp_body(:,2);
z_p = Xp_body(:,3);
h_nplus1 = plot3(x_p, y_p, z_p, 'w', 'LineWidth', 2);
plot3(x_p(1), y_p(1), z_p(1), 'wo', 'MarkerFaceColor', 'w');
plot3(x_p(end), y_p(end), z_p(end), 'ws', 'MarkerFaceColor', 'w');

for i = 1:n
    
    x_i = X(:,3*(i-1)+1);
    y_i = X(:,3*(i-1)+2);
    z_i = X(:,3*(i-1)+3);
    
    % Trajectory
    plot3(x_i, y_i, z_i, 'Color', colors(i,:), 'LineWidth', 1.5);
    
    % Initial position
    plot3(x_i(1), y_i(1), z_i(1), 'o', ...
        'Color', colors(i,:), 'MarkerFaceColor', colors(i,:));
    
    % Final position
    plot3(x_i(end), y_i(end), z_i(end), 's', ...
        'Color', colors(i,:), 'MarkerFaceColor', colors(i,:));
    
    % Body label (clean + readable)
    text(x_i(end), y_i(end), z_i(end), ...
        ['  ' names{i}], ...
        'Color', colors(i,:), ...
        'FontSize', 10, ...
        'FontWeight', 'bold');
end

grid on
axis equal
xlabel('x (AU)')
ylabel('y (AU)')
zlabel('z (AU)')
title('3D N+1st-Body Bodycentric Trajectories (AU, Years, Solar Masses)')

legend([h_nplus1, h_init, h_final],{'N+1st Body -- Satellite','Initial Position','Final Position'},'location','southeast')

view(3)

Nt = length(tq_body);

% Build interpolants
r_interp_sc = griddedInterpolant(tq_bary, ...
    Xp_bary(:,1:3), ...
    'linear', 'nearest');

v_interp_sc = griddedInterpolant(tq_bary, ...
    Xp_bary(:,4:6), ...
    'linear', 'nearest');

r_interp_N = griddedInterpolant(t_NBP_bary, ...
    X_NBP_bary(:,3*(idxN-1)+1:3*idxN), ...
    'linear', 'nearest');

v_interp_N = griddedInterpolant(t_NBP_bary, ...
    X_NBP_bary(:,3*n + 3*(idxN-1)+1:3*n + 3*idxN), ...
    'linear', 'nearest');

err = zeros(Nt,1);
err_v = zeros(Nt,1);

for k = 1:Nt

    % Bodycentric simulation solution
    r_p = Xp_body(k,1:3)';
    v_p = Xp_body(k,4:6)';

    % Interpolate barycentric solution
    r_p_bary = r_interp_sc(tq_body(k))';
    v_p_bary = v_interp_sc(tq_body(k))';

    % Interpolate central body from N-body
    r_N = r_interp_N(tq_body(k))';
    v_N = v_interp_N(tq_body(k))';

    % Convert barycentric (truth) → bodycentric
    r_true_p = r_p_bary - r_N;
    v_true_p = v_p_bary - v_N;

    % Compute error in same frame
    err(k)   = norm(r_true_p - r_p);
    err_v(k) = norm(v_true_p - v_p);

end
% Report max error
fprintf('Max relative error between bodycentric and barycentric position (N+1 restricted): %.3e AU\n', max(err));
fprintf('Max relative error between bodycentric and barycentric velocity (N+1 restricted): %.3e AU/day\n', max(err_v));

% The restricted (N+1)-body bodycentric simulation agrees closely with the 
% barycentric formulation. The maximum position and velocity differences 
% between the two approaches are on the order of 2.4-5 AU and 7.5e-5 AU/day, 
% respectively. This demonstrates that both formulations are dynamically 
% equivalent and that the bodycentric implementation correctly reproduces 
% the particle motion relative to the central body.