function dXdt = restricted_NBP_barycentric(tq, Xp, mass_bodies, t, X, G)

% tq = current time (from ode45)
% Xp = [r_p; v_p] for particle

r_p = Xp(1:3);
v_p = Xp(4:6);

n = length(mass_bodies);
a_p = zeros(3,1);

for j = 1:n
    
    m_j = mass_bodies(j);

    % Extract trajectory of body j
    rj_hist = X(:, 3*(j-1)+1 : 3*j);   % (Nt x 3)

    % Interpolate position at current time tq
    r_j = interp1(t, rj_hist, tq, 'linear', 'extrap')';    % (3x1)

    % Relative vector
    r_pj = r_j - r_p;

    % Acceleration contribution
    a_p = a_p + (G*m_j / norm(r_pj)^3) * r_pj;
end

dXdt = [v_p; a_p];

end