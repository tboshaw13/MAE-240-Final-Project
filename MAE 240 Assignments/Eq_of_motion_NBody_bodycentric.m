function dXdt = Eq_of_motion_NBody_bodycentric(~, X, mass_bodies, G)
% INPUTS
% t -- time vector 1xlength(t)
% X -- Trajectory vector (pos,vel of all bodies)
% mass_bodies -- mass of all bodies
% G -- gravitational constant

% OUTPUTS
% dXdt -- change of trajectory over time step for numerical integration


n = numel(mass_bodies);       % # of bodies in system

r = X(1:3*n);     % Position of body trajectory
v = X(3*n+1:end); % Velocity of body trajectory

% Make motion relative to largest body
[m_N,idxN] = max(mass_bodies); % Largest body 'Nth" body
r_N = r(3*(idxN-1)+1:3*idxN); % Position of "Nth" body's trajectory

% Initialize the output derivative vector
dvdt = zeros(3*n,1);

% Compute accelerations due to gravitational forces
for i = 1:n % ensures looping through all bodies

    if i ~= idxN

        r_i = r(3*(i-1)+1:3*i); % Position of "i" body's trajectory
        m_i = mass_bodies(i);   % mass of ith body
        script_r_i = r_i - r_N; % redefine postion from Nth body

        dvdt_i = - ( G*(m_N+m_i)*script_r_i*norm(script_r_i)^-3 ); % Direct term
        
        for j = 1:n   % ensures looping through all bodies
            if i ~= j && j ~= idxN % ignore whenever j ==  && j == N
    
                m_j = mass_bodies(j);          % mass of "j" body
                r_j = r(3*(j-1)+1:3*j);        % Position of "j" body's trajectory
                script_r_j = r_j - r_N;        % Position of jth body wrt Nth body
                script_r_ji = r_i - r_j;       % Difference in body's relative position
                
                % Summation of all other bodies' grav. pulls on "ith" body
                dvdt_i = dvdt_i - ( (G * m_j) * ((norm(script_r_ji)^-3 * script_r_ji) + (norm(script_r_j)^-3 * script_r_j)) ); % Indirect term
            end
        end

        % Compute Derivative vector
        dvdt(3*(i-1)+1:3*i) = dvdt_i; % acceleration of "ith" body
    end    
end

% Derivative of state vectors [dr/dt; dv/dt]
dXdt = [v; dvdt];

end