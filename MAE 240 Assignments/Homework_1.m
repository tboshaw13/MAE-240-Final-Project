%% Problem 4 -- N-Body-Problem 
% Part a) N-Body Solar System Simulation (AU, years, solar masses)

clear; clc; close all;

% Constants
G = 4*pi^2;   % AU^3 / (yr^2 * solar mass)

% I chose the Sun, Eath, Mars, and Jupiter as my masses becuase their
% masses, orbital period, and relative positions to one another are known.
% This means my simulation has some knowledge of what should happen,
% allowing me to check my work. The inital positions all started on the x
% axis, with the time period being one year (a full rotation around the sun
% for Earth). The initial velocities were set as the units over time, which
% means that since the Earth travels 2*pi AU per year, everything else can
% be normalized by these units. 

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

% Combine into state vector
X0 = [r0; v0];

% Time span (years)
tspan = [0 20];

% ODE Solver options for stability
options = odeset('RelTol',1e-10,'AbsTol',1e-12);

% Run N-body problem eq. of motion through solver
[t, X] = ode45(@(t,X) Eq_of_motion_NBody(t,X,mass_bodies,G), tspan, X0, options);
n = numel(mass_bodies);

% 3D Trajectory Plot
names = {'Sun','Earth','Mars','Jupiter'};
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
title('3D N-Body Trajectories (AU, Years, Solar Masses)')

legend([h_init, h_final], ...
    {'Initial Position','Final Position'})

view(3)

% As seen in the graph below, I chose to represent 4 Bodies in
% the solar system. The Sun, the Earth, Mars, and Jupiter. The initial
% positions are labeled as circles and the final positions as squares.

%% Part b) Verify Problem 1 -- COM Prediction accuracy

% Center of Mass Calculation
M = sum(mass_bodies); % total mass
r_c = zeros(length(t),3);

for k = 1:length(t)
    
    r_c_k = zeros(3,1);
    
    for i = 1:n
        r_i = X(k,3*(i-1)+1 : 3*i)'; % position
        r_c_k = r_c_k + mass_bodies(i)*r_i;
    end
    
    r_c(k,:) = (r_c_k / M)';
end

% Initial COM Velocity (at t0)
v_c0 = zeros(3,1);

for i = 1:n
    v_i = X(1,3*n + 3*(i-1)+1 : 3*n + 3*i)'; % velocity
    v_c0 = v_c0 + mass_bodies(i)*v_i;
end

v_c0 = v_c0 / M; % divided by sum of all masses

% Euler: Linear prediction
r_c_pred = zeros(size(r_c));

r_c0 = r_c(1,:)'; % rc0 = rc(t0)

for k = 1:length(t)
    r_c_pred(k,:) = (r_c0 + v_c0*(t(k) - t(1)))'; % rc(t) = rc(0)+vc(t-t0)
end

% Plot COM comparison (difference)

figure; hold on;

plot3(r_c(:,1), r_c(:,2), r_c(:,3), 'b', 'LineWidth', 2)
plot3(r_c_pred(:,1), r_c_pred(:,2), r_c_pred(:,3), 'r--', 'LineWidth', 2)

% Initial point
plot3(r_c(1,1), r_c(1,2), r_c(1,3), 'go', ...
    'MarkerFaceColor','g', 'MarkerSize',8)

% Final point
plot3(r_c(end,1), r_c(end,2), r_c(end,3), 'ws', ...
    'MarkerFaceColor','w', 'MarkerSize',16)

% Zoom in 
axis equal
xlim([r_c(1,1)-1e-15, r_c(1,1)+1e-15])
ylim([r_c(1,2)-1e-15, r_c(1,2)+1e-15])
zlim([r_c(1,3)-1e-15, r_c(1,3)+1e-15])


grid on

xlabel('x (AU)')
ylabel('y (AU)')
zlabel('z (AU)')
title('Center of Mass Motion')

legend('Numerical COM','Linear Prediction','Initial','Final')
view(3)

% Error calculation
error = vecnorm(r_c - r_c_pred, 2, 2); % Euclidean norm at each time
max_error = max(error); % find maxima

disp(['Max COM deviation: ', num2str(max_error)])

% The graph below shows the COM as calculated numerically, as well as with
% the linear prediction. It is easy to see that the COM barely moves (as
% the sun's mass dominates the other planets), but when it does, the linear
% prediction follows the same line. This shows that the linear prediction
% is very accurate for this example. 

%% Part c) Verify Problem 2 -- Angular Momentum Calculation

% Preallocate
H = zeros(length(t),3);

% Loop through calculated position and velocities of the bodies
for k = 1:length(t)
    
    % Preallcoate H
    H_k = zeros(3,1);
    
    for i = 1:n
        
        r_i = X(k,3*(i-1)+1 : 3*i)'; % position at every "ith" body
        v_i = X(k,3*n + 3*(i-1)+1 : 3*n + 3*i)'; % velocity at every "ith" body
        
        H_k = H_k + mass_bodies(i) * cross(r_i, v_i); % H = mi * (ri x vi)
    end
    
    H(k,:) = H_k';
end

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

% Relative Drift
H0 = H(1,:); % initial angular momentum 

% norm of the difference between the inital and actual angular momentum
delta_H = vecnorm(H - H0, 2, 2); 
rel_drift = delta_H / norm(H0); % see how much angular momentum is changing

max_rel_drift = max(rel_drift); % find maximum drift

disp(['Max relative angular momentum drift: ', num2str(max_rel_drift)])

% Planarity Check

z_all = X(:,3:3:3*n); % all z components of pos. & vel.

max_z = max(abs(z_all), [], 'all'); % find maxima

disp(['Max deviation from plane (z): ', num2str(max_z)])

% The plot below shows each component (x,y,z) of the angular momentum of
% the 4 Body system over time. It is easy to see that the x and y
% components are constant at 0, meaning the system is a planar system with
% the direction of angular momentum (here along the z axis) is
% perpendicular to the motion of the bodies. The z component is shown to be
% constant, with a very small drift with a magnitude of e-13. The planarity
% check was a numerical check that shows that every z component of position
% and velocity was 0, which means that every body's position and velocity
% was on the x-y plane.

%% Part d) Verify Problem 3 -- Energy Calculation

% Preallocate Kinetic and Potential Energy 
T = zeros(length(t),1);
U = zeros(length(t),1);

for k = 1:length(t)
    
    T_k = 0;
    U_k = 0;
    
    % Calculate Kinetic Energy 
    for i = 1:n
        v_i = X(k,3*n + 3*(i-1)+1 : 3*n + 3*i); % velocity 
        T_k = T_k + 0.5 * mass_bodies(i) * dot(v_i, v_i); % 0.5*mv^2
    end
    
    % Calculate Potential Energy 
    for i = 1:n
        r_i = X(k,3*(i-1)+1 : 3*i); % position of "ith" body
        
        for j = i+1:n
            r_j = X(k,3*(j-1)+1 : 3*j); % position of the other 3 bodies
            r_ij = norm(r_i - r_j); % difference in position
            
            % U = G*mi*mj/rij
            U_k = U_k + G * mass_bodies(i)*mass_bodies(j) / r_ij; 
        end
    end
    
    T(k) = T_k;
    U(k) = U_k;
end

E = T - U; % Total energy

% Normalized Energy Plot

E0 = abs(E(1)); % compute normalized en. at inital time t0

figure; hold on;

plot(t, T/E0, 'r', 'LineWidth', 1.5) % Kinetic Energy plot
plot(t, U/E0, 'g', 'LineWidth', 1.5) % Potential Energy plot
plot(t, E/E0, 'b', 'LineWidth', 2) % Total Energy plot

xlabel('Time (years)')
ylabel('Normalized Energy')
title('Energy Components (Normalized by |E(t0)|)')

legend('T/|E_0|','U/|E_0|','E/|E_0|')
grid on

% Relative Energy Error plot

% Calculate relative error in energy over time
rel_error = abs(E - E(1)) / abs(E(1));

figure;
semilogy(t, rel_error, 'LineWidth', 2)

xlabel('Time (years)')
ylabel('|E(t) - E(0)| / |E(0)|')
title('Relative Energy Error (Log Scale)')
grid on

% The first graph below shows the potential, kinetic, and total energy normalized
% by the total energy. It can be easily seen that the kinetic, potential,
% and total energies are all constant over time, with the potential energy
% being 2x the kinetic energy. 

% As seen in the second graph below, the energy error (or deviation from constant
% energy over time) varies from e-16 to e-11, which is basically just
% floating point error. This shows that the energy is accurately constant
% over time. 