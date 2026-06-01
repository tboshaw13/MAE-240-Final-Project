function [a, e, i, cap_omega, omega, M] = Cartesian_to_oe(X0, mu, delta_t)

% Extract position and velocity
r_vec = X0(1:3);
v_vec = X0(4:6);

r = norm(r_vec);
v = norm(v_vec);

r_hat = r_vec/r;

% Specific orbital energy
E = v^2/2 - mu/r;

% Angular momentum
H_vec = cross(r_vec,v_vec);
H = norm(H_vec);

if H > 1e-12
    H_hat = H_vec/H;
else
    H_hat = [NaN;NaN;NaN];
end

% Eccentricity vector
e_vec = (1/mu)*cross(v_vec,H_vec) - r_hat;
e = norm(e_vec);

if e > 1e-12
    e_hat = e_vec/e;
else
    e_hat = [NaN;NaN;NaN];
end

% Reference axes
x_hat = [1;0;0];
z_hat = [0;0;1];

% Node vector
n_vec = cross(z_hat,H_vec);
n = norm(n_vec);

if n > 1e-12
    n_hat = n_vec/n;
else
    n_hat = [NaN;NaN;NaN];
end

% Semi-major axis
if abs(E) > 1e-12
    a = -mu/(2*E);
else
    a = Inf;
end

% Inclination
if H > 1e-12
    arg = dot(H_hat,z_hat);
    arg = max(-1,min(1,arg));
    i = acos(arg);
else
    i = NaN;
end

% RAAN (cap_omega)
if n > 1e-12

    arg = max(-1,min(1,n_hat(1)));

    cap_omega = acos(arg);

    if n_vec(2) < 0
        cap_omega = 2*pi-cap_omega;
    end

else

    cap_omega = NaN;

end

% Argument of periapsis
if e > 1e-12 && n > 1e-12

    omega = atan2( ...
        dot(H_hat,cross(n_hat,e_vec)), ...
        dot(n_hat,e_vec));

    omega = mod(omega,2*pi);

elseif e > 1e-12 && n <= 1e-12

    arg = max(-1,min(1,e_hat(1)));

    omega = acos(arg);

    if e_vec(2) < 0
        omega = 2*pi-omega;
    end

else

    omega = NaN;

end

% True anomaly
if e > 1e-12

    arg = dot(r_hat,e_hat);
    arg = max(-1,min(1,arg));

    f = acos(arg);

    if dot(r_vec,v_vec) < 0
        f = 2*pi-f;
    end

else

    f = NaN;

end

% Eccentric anomaly
if e < 1

    E_anom = atan2( ...
        sqrt(max(0,1-e^2))*sin(f), ...
        e+cos(f));

else

    E_anom = NaN;

end

% Mean motion
if a > 0 && isfinite(a)

    n_mean = sqrt(mu/a^3);

else

    n_mean = NaN;

end

% Mean anomaly
if ~isnan(E_anom)

    M0 = E_anom - e*sin(E_anom);

    M = M0 + n_mean*delta_t;

    M = mod(M,2*pi);

else

    M = NaN;

end

% Convert to degrees
i = rad2deg(i);
cap_omega = rad2deg(cap_omega);
omega = rad2deg(omega);
M = rad2deg(M);

end