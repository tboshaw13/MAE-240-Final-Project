function dXdt = restricted_NBP_bodycentric(tq, Xp, mass_bodies, r_interp, G)
% INPUTS
% tq -- current time
% Xp -- [r_p; v_p] particle (bodycentric)
% mass_bodies -- masses
% r_interp -- interpolants of body trajectories
% G -- gravitational constant

r_p = Xp(1:3);
v_p = Xp(4:6);

n = length(mass_bodies);

[m_N, idxN] = max(mass_bodies);

% Particle already expressed in bodycentric frame
script_r_p = r_p;

% Direct term
dvdt_p = -G*m_N * script_r_p / max(norm(script_r_p)^3, eps);

for j = 1:n
    if j ~= idxN
    
        m_j = mass_bodies(j);
        
        % Interpolate barycentric position of body j
        r_j = r_interp{j}(tq)';
        
        % Interpolate barycentric position of central body
        r_N = r_interp{idxN}(tq)';
        
        % Convert body j to bodycentric frame
        script_r_j = r_j - r_N;
        
        % Relative particle vector
        script_r_jp = script_r_p - script_r_j;
        
        % Indirect terms
        dvdt_p = dvdt_p - G*m_j * ( ...
            script_r_jp / max(norm(script_r_jp)^3, eps) + ...
            script_r_j / max(norm(script_r_j)^3, eps) );
    end
end

dXdt = [v_p; dvdt_p];

end