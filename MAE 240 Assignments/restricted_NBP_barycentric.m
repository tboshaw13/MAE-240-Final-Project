function dXdt = restricted_NBP_barycentric(tq, Xp, mass_bodies, r_interp, G)

% tq = current time (from ode45)
% Xp = [r_p; v_p] for particle

r_p = Xp(1:3);
v_p = Xp(4:6);

n = length(mass_bodies);
dvdt_p = zeros(3,1);

for j = 1:n
    
    m_j = mass_bodies(j);

    % Interpolate position of body j at current time tq
    r_j = r_interp{j}(tq)';   % (3x1)

    % Relative vector
    r_jp = r_p - r_j;

    % Acceleration contribution
    dvdt_p = dvdt_p - G*m_j * r_jp / max(norm(r_jp)^3, eps);

end

dXdt = [v_p; dvdt_p];

end