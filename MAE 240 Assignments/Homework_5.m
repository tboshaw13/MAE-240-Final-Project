%% Homework 4: Question 3
clc;clear;clf;

% Define Initial Conditions
delta_t = 0;
mu = 398600.4418;  % Gravitational parameter for problem (km^3/s^2)

% Orbital Elements for 3 Orbits

a    = [7000, 12000, 26560];      % km
e    = [0.10, 0.25, 0.15];
i          = [45, 70, 55];        % deg
cap_omega  = [60, 120, 250];      % deg
omega      = [30, 75, 110];       % deg
M0         = [0, 45, 180];        % deg
ad = [0; 0; 1e-5];               % km/s Acceleration due to other forces (if any)

for k = 1:length(a)
    T = 2*pi*sqrt(a(k)^3/mu);
    tspan = [0 10*T];
    
    [X0] = oe_to_cartesian(a(k), e(k), i(k), cap_omega(k), omega(k), M0(k), mu, delta_t);
    
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
    
    % Compute resultant position and velocity magnitudes
    r_mag = vecnorm(X_out(:,1:3),2,2);
    v_mag = vecnorm(X_out(:,4:6),2,2);
    
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
    for j = 1:N
    
        Xj = X_out(j,:)';
    
        [a_hist(j), ...
         e_hist(j), ...
         i_hist(j), ...
         cap_omega_hist(j), ...
         omega_hist(j), ...
         M_hist(j)] = ...
            Cartesian_to_oe(Xj,mu,t_out(j));
    
    end

    % Theory Part c)
    
    n = sqrt(mu/a(k)^3);
    E_hist = zeros(N,1);
    
    a_theory = zeros(N,1);
    e_theory = zeros(N,1);
    
    for j = 1:N
    
        % Unperturbed mean anomaly evolution
        Mj = M0(k) + n*t_out(j);
    
        % Wrap angle
        Mj = mod(Mj,2*pi);
    
        % Solve Kepler's Equation
        keplerEq = @(E) E - e(k)*sin(E) - Mj;
    
        Ej = fzero(keplerEq,Mj);
        E_hist(j) = Ej;
    
        % Delta-a
        delta_a = (2*sin(i(k))/n^2) * ...
            ( sqrt(1 - e(k)^2)*cos(omega(k))*sin(Ej) ...
            - sin(omega(k))*(1 - cos(Ej)) ) ...
            * ad(3);
    
        % Delta-e
        delta_e = (sqrt(1 - e(k)^2)*sin(i(k))/(n^2*a(k))) * ...
            ( (3/2)*cos(omega(k))*Ej ...
            - 2*e(k)*cos(omega(k))*sin(Ej) ...
            + (1/4)*cos(omega(k))*sin(2*Ej) ...
            - (sqrt(1 - e(k)^2)/4)*sin(omega(k))*(1 - cos(2*Ej)) ) ...
            * ad(3);
    
        a_theory(j) = a(k) + delta_a;
        e_theory(j) = e(k) + delta_e;
    
    end
    
    % Plot Theory vs Numerical Separately
    
    figure(k)
    
    tiledlayout(2,2,'TileSpacing','compact')
    
    % Theory Semi-Major Axis
    nexttile
    
    plot(t_out,a_theory,'LineWidth',1.5)
    
    grid on
    
    ylabel('a [km]')
    
    title('Theory Semi-Major Axis')
    
    % Numerical Semi-Major Axis
    nexttile
    
    plot(t_out,a_hist,'LineWidth',1.5)
    
    grid on
    
    ylabel('a [km]')
    
    title('Numerical Semi-Major Axis')
    
    % Theory Eccentricity
    nexttile
    
    plot(t_out,e_theory,'LineWidth',1.5)
    
    grid on
    
    ylabel('e')
    xlabel('Time [s]')
    
    title('Theory Eccentricity')
    
    % Numerical Eccentricity
    nexttile
    
    plot(t_out,e_hist,'LineWidth',1.5)
    
    grid on
    
    ylabel('e')
    xlabel('Time [s]')
    
    title('Numerical Eccentricity')
    
    sgtitle(['Orbit ',num2str(k),': Theory vs Numerical'])

    % Theory part d)
    % Mean Short-Period Values

    % Mean of Delta-a short-period term
    mean_delta_a_theory = mean(a_theory);
    
    % Remove secular eccentricity term
    delta_e_sec = (sqrt(1 - e(k)^2)*sind(i(k))/(n^2*a(k))) .* ...
        ((3/2)*cosd(omega(k)).*E_hist) * ad(3);
    
    % Short-period eccentricity only
    delta_e_short = e_theory - delta_e_sec;
    
    % Mean short-period eccentricity
    mean_delta_e_theory = mean(delta_e_short);
    
    % Numerical Means
    
    % Numerical mean semi-major axis shift
    mean_delta_a_num = mean(a_hist);
    
    % Numerical mean eccentricity shift
    mean_delta_e_num = mean(e_hist);
    
    % Display Results
    
    fprintf('\nOrbit %d\n',k)
    
    fprintf('Theory Mean Delta-a      = %.6e km\n',mean_delta_a_theory)
    fprintf('Numerical Mean Delta-a  = %.6e km\n',mean_delta_a_num)
    
    fprintf('Theory Mean Delta-e(sp) = %.6e\n',mean_delta_e_theory)
    fprintf('Numerical Mean Delta-e  = %.6e\n',mean_delta_e_num)
end

