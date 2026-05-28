function [dXdt,R] = Eq_of_motion_NBody_Barycentric(~, X, mass_bodies, G)
% INPUTS
% t -- time vector 1xlength(t)
% X -- Trajectory vector (pos,vel of all bodies)
% mass_bodies -- mass of all bodies
% G -- gravitational constant

% OUTPUTS
% dXdt -- change of trajectory over time step for numerical integration


n = numel(mass_bodies); % # of bodies in system
M = sum(mass_bodies);   % total mass of all bodies

r = X(1:3*n);     % Position of body trajectory
v = X(3*n+1:end); % Velocity of body trajectory
R = 0;

% Initialize the output derivative vector
dvdt = zeros(3*n,1);
for i = 1:n % ensures looping through all bodies
    r_i = r(3*(i-1)+1:3*i); % Position of "i" body's trajectory
    R = mass_bodies(i)*r_i + R;
end
R = (1/M) * R; % Barycenter (COM)

% Compute accelerations due to gravitational forces
for i = 1:n % ensures looping through all bodies

    r_i = r(3*(i-1)+1:3*i); % Position of "i" body's trajectory
    dvdt_i = zeros(3,1);    % Preallocate derivative vector for every "i" body

    for j = 1:n   % ensures looping through all bodies
        if i ~= j % ignore whenever i == j

            m_j = mass_bodies(j);   % mass of "j" body
            r_j = r(3*(j-1)+1:3*j); % Position of "j" body's trajectory
            % r_i_script = r_i - R;   % position of r_i relative to barycenter
            % r_j_script = r_j - R;   % position of r_j relative to barycenter
            r_ij = r_j - r_i;       % Difference in body's 
            
            % Summation of all other bodies' grav. pulls on "ith" body
            dvdt_i = dvdt_i + (((G * m_j) / norm(r_ij)^3) * r_ij);
        end
    end
    % Compute Derivative vector
    dvdt(3*(i-1)+1:3*i) = dvdt_i; % acceleration of "ith" body
end


% Derivative of state vectors [dr/dt; dv/dt]
dXdt = [v; dvdt];
end