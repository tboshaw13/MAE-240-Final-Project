% Dynamics Function
function Xdot = twoBodyConstAccelODE(~, X, mu, ad)

    % Extract position and velocity
    r = X(1:3);
    v = X(4:6);

    % Radius magnitude
    r_norm = norm(r);

    % Two-body acceleration
    a_grav = -mu/r_norm^3 * r;

    % Total acceleration
    a_total = a_grav + ad(:);

    % State derivative
    Xdot = [v;
            a_total];

end