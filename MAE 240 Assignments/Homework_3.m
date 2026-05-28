%% Homework 3 -- Cartesian -> OE Test Problem 1
clc;clear;

% Define parameters for the Cartesian to OE transformation
mu = 398600.4354; % Gravitational parameter for Earth in km^3/s^2
X0 = [-1264.608 8013.809 -3371.252 -6.039621 -0.204398 2.096715]'; % Trajectory [r;v] ~ [km;km/s]
delta_t = 3600; % [s]

% Call Cartesian -> OE function
[a, e, i, cap_omega, omega, M] = Cartesian_to_oe(X0, mu, delta_t);

fprintf('Orbital Elements:\n');
fprintf('a         = %.6f\n', a);
fprintf('e         = %.6f\n', e);
fprintf('i         = %.6f deg\n', i);
fprintf('Omega     = %.6f deg\n', cap_omega);
fprintf('omega     = %.6f deg\n', omega);
fprintf('M         = %.6f deg\n\n', M);

% Inital M0
M0 = 10;

% Call OE -> Cartesian function
[X0] = oe_to_cartesian(a, e, i, cap_omega, omega, M0, mu, delta_t);
fprintf('R         = %.6f km\n', X0(1:3));
fprintf('V         = %.6f km/s\n', X0(4:6));
