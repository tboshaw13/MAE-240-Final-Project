%% Homework 4: Question 3
clc;clear;clf;

% Define Initial Conditions
% tspan = [0, 60*24*60*50]; % Time span for propagation [s]
delta_t = 0;
mu = 398600.4418;  % Gravitational parameter for problem (km^3/s^2)

% Orbit 1
a     = 7000;        % km
e     = 0.10;
i     = (45);        % rad
cap_omega  = (60);   % rad
omega     = (30);    % rad
M0    = (0);         % rad
T1 = 2*pi*sqrt(a^3/mu);
% T2 = 2*pi*sqrt(a2^3/mu);
% T3 = 2*pi*sqrt(a3^3/mu);

tspan = [0 20*T1];
% tspan2 = [0 20*T2];
% tspan3 = [0 20*T3];
[X0] = oe_to_cartesian(a, e, i, cap_omega, omega, M0, mu, delta_t);
r0 = X0(1:3);
v0 = X0(4:6);

% r0 = [7000;0;0];   % Initial position vector [km] 
% v0 = [0;5.5;5.5];  % Initial velocity vector [km/s]
% tspan = [0, 20*pi]; % Time span for propagation [s]
% delta_t = 0;
% mu = 1;  % Gravitational parameter for problem (km^3/s^2)
% r0 = [1;0;0];   % Initial position vector [km] 
% v0 = [0;1;0];  % Initial velocity vector [km/s]
%a = [0.001, 0.01, 0.1, 1];
a = [10e-7, 0.01, 0.1, 1];

% for i = 1:length(a)
for i = 1
    ad = [0; 0; a(i)];    % Acceleration due to other forces (if any)
    
    % Initial state vector
    X0 = [r0(:); v0(:)];
    
    % Integrate with ode45
    options = odeset('RelTol',1e-12,'AbsTol',1e-12);
    
    [t_out, X_out] = ode45(@(t,X) twoBodyConstAccelODE(t,X,mu,ad),tspan, X0, options);
    
    % Extract states
    rx  = X_out(:,1); % Position
    ry  = X_out(:,2);
    rz  = X_out(:,3);
    vx = X_out(:,4);  % Velocity
    vy = X_out(:,5);
    vz = X_out(:,6);
    
    % % Plot Position and Velocity as a funciton of time
    % figure;
    % plot3(rx,ry,rz,'LineWidth',1.5);
    % grid on;
    % axis equal;
    % xlabel('x [km]');
    % ylabel('y [km]');
    % zlabel('z [km]');
    % title('3D Orbit Trajectory');
    
    % Compute resultant position and velocity magnitudes
    r_mag = vecnorm(X_out(:,1:3),2,2);
    v_mag = vecnorm(X_out(:,4:6),2,2);
    
    % % Plot resultant position and velocity vs time
    % figure; 
    % plot(t_out,r_mag,'LineWidth',1.5)
    % ylabel('|r| [km]')
    % xlabel('Time [s]')
    % title('Position vs Time')
    % 
    % figure;
    % plot(t_out,v_mag,'LineWidth',1.5)
    % ylabel('|v| [km/s]')
    % xlabel('Time [s]')
    % title('Velocity vs Time')
    % grid on
    
    %Compute Classical Orbital Elements
    % Preallocate orbital element histories
    N = length(t_out);
    
    a_hist         = zeros(N,1);
    e_hist         = zeros(N,1);
    i_hist         = zeros(N,1);
    cap_omega_hist = zeros(N,1);
    omega_hist     = zeros(N,1);
    M_hist         = zeros(N,1);
    
    % Compute OE history
    for k = 1:N
    
        Xk = X_out(k,:)';
    
        [a_hist(k), ...
         e_hist(k), ...
         i_hist(k), ...
         cap_omega_hist(k), ...
         omega_hist(k), ...
         M_hist(k)] = ...
            Cartesian_to_oe(Xk,mu,t_out(k));
    
    end
    
    % % Plot OE over time
    % figure
    % 
    % % Semi-major axis
    % subplot(3,2,1)
    % plot(t_out,a_hist,'LineWidth',1.5)
    % grid on
    % ylabel('a [km]')
    % title('Semi-Major Axis')
    % 
    % % Eccentricity
    % subplot(3,2,2)
    % plot(t_out,e_hist,'LineWidth',1.5)
    % grid on
    % ylabel('e')
    % title('Eccentricity')
    % 
    % % Inclination
    % subplot(3,2,3)
    % plot(t_out,i_hist,'LineWidth',1.5)
    % grid on
    % ylabel('i [deg]')
    % title('Inclination')
    % 
    % % % RAAN
    % % subplot(3,2,4)
    % % plot(t_out,cap_omega_hist,'LineWidth',1.5)
    % % grid on
    % % ylabel('\Omega [deg]')
    % % title('RAAN')
    % 
    % % Argument of Periapsis
    % subplot(3,2,5)
    % plot(t_out,omega_hist,'LineWidth',1.5)
    % grid on
    % ylabel('\omega [deg]')
    % xlabel('Time [s]')
    % title('Argument of Periapsis')
    % 
    % % Mean Anomaly
    % subplot(3,2,6)
    % plot(t_out,M_hist,'LineWidth',1.5)
    % grid on
    % ylabel('M [deg]')
    % xlabel('Time [s]')
    % title('Mean Anomaly')
    % 
    % sgtitle('Classical Orbital Elements vs Time')

    figure (i)
    
    tiledlayout(3,4,'TileSpacing','compact')
    
    % Plot Orbit Trajectory
    nexttile([3 2])
    
    plot3(rx,ry,rz,'LineWidth',1.5)
    grid on  
    xlabel('x')
    ylabel('y')
    zlabel('z')
    
    title(['3D Orbit Trajectory, a_d = ',num2str(a(i))])
    
    % Plot Semi-Major Axis
    nexttile
    
    plot(t_out,a_hist,'LineWidth',1.2)
    
    grid on
    
    ylabel('a')
    
    title('Semi-Major Axis')
    
    % Plot Eccentricity
    nexttile
    
    plot(t_out,e_hist,'LineWidth',1.2)
    
    grid on
    
    ylabel('e')
    
    title('Eccentricity')
    
    % Plot Inclination
    nexttile
    
    plot(t_out,i_hist,'LineWidth',1.2)
    
    grid on
    
    ylabel('i [deg]')
    
    title('Inclination')
    
    % Plot Argument of Periapsis
    nexttile
    
    plot(t_out,omega_hist,'LineWidth',1.2)
    
    grid on
    
    ylabel('\omega [deg]')
    
    title('Argument of Periapsis')
    
    % Plot Mean Anomaly
    nexttile
    
    plot(t_out,M_hist,'LineWidth',1.2)
    
    grid on
    
    ylabel('M [deg]')
    
    xlabel('Time')
    
    title('Mean Anomaly')
    
    % Plot Empty Tile
    nexttile
    axis off
    
    % Overall Figure Title
    sgtitle(['Constant Disturbance Acceleration: a_d = ',num2str(a(i))])
end

