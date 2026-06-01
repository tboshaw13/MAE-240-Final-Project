function [X0] = oe_to_cartesian(a, e, i, cap_omega, omega, M0, mu, delta_t)

% Convert degrees → radians
i = deg2rad(i);
cap_omega = deg2rad(cap_omega);
omega = deg2rad(omega);
M0 = deg2rad(M0);

% Mean motion
n = sqrt(mu / a^3);

% Propagate mean anomaly
M = M0 + n * delta_t;
M = mod(M, 2*pi);

% Solve Kepler's equation: M = E - e sin(E)
E = M; % initial guess
for k = 1:10
    E = M + e*sin(E);
end

% True anomaly
f = 2 * atan2( sqrt(1+e)*sin(E/2), sqrt(1-e)*cos(E/2) );

% Radius
r_mag = a * (1 - e*cos(E));

% Perifocal position
r_perifocal = [r_mag*cos(f);
               r_mag*sin(f);
               0];

% Perifocal velocity
v_perifocal = sqrt(mu*a)/r_mag * ...
              [-sin(E);
                sqrt(1-e^2)*cos(E);
                0];

% Rotation matrix (PQW → ECI)
R = [cos(cap_omega)*cos(omega) - sin(cap_omega)*sin(omega)*cos(i), ...
    -cos(cap_omega)*sin(omega) - sin(cap_omega)*cos(omega)*cos(i), ...
     sin(cap_omega)*sin(i);

     sin(cap_omega)*cos(omega) + cos(cap_omega)*sin(omega)*cos(i), ...
    -sin(cap_omega)*sin(omega) + cos(cap_omega)*cos(omega)*cos(i), ...
    -cos(cap_omega)*sin(i);

     sin(omega)*sin(i), ...
     cos(omega)*sin(i), ...
     cos(i)];

% Convert to inertial frame
r = R * r_perifocal;
v = R * v_perifocal;

X0 = [r;v];

end