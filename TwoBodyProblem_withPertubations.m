function dydt = TwoBodyProblem_withPertubations(t, X, const)

    % Extract position and velocity from state vector
    rVec = X(1:3);
    vVec = X(4:6);

    x = rVec(1);
    y = rVec(2);
    z = rVec(3);

    vx = vVec(1);
    vy = vVec(2);
    vz = vVec(3);

    % Magnitudes
    r = norm(rVec);
    v = norm(vVec);

    %% Two-body Earth gravity
    a_grav = -const.mu * rVec / r^3;

    %% J2 perturbation
    J2_factor = -3 * const.mu * const.J2 * const.R^2 / (2 * r^4);

    ax_J2 = J2_factor * (x/r) * (1 - 5 * (z^2 / r^2));
    ay_J2 = J2_factor * (y/r) * (1 - 5 * (z^2 / r^2));
    az_J2 = J2_factor * (z/r) * (3 - 5 * (z^2 / r^2));

    a_J2 = [ax_J2; ay_J2; az_J2];

    % %% Atmospheric drag
    % % For GEO, this is basically negligible, but keep it if you want to compare models.
    % h = r - const.R; % altitude, km
    % [rho, ~] = atmosphere(h); % rho should be kg/m^3
    % 
    % % v is in km/s, so convert velocity to m/s for drag, then back to km/s^2
    % v_mps = v * 1000;
    % vVec_mps = vVec * 1000;
    % 
    % a_drag_mps2 = -0.5 * rho * const.Cd * const.A_m * v_mps * vVec_mps;
    % a_drag = a_drag_mps2 / 1000; % convert m/s^2 to km/s^2

    %% Approximate Sun and Moon position vectors relative to Earth
    rSun = sunPositionApprox(t, const);    % km
    rMoon = moonPositionApprox(t, const); % km

    %% Third-body gravity from Sun and Moon
    a_Sun_grav = thirdBodyAccel(rVec, rSun, const.muSun);
    a_Moon_grav = thirdBodyAccel(rVec, rMoon, const.muMoon);

    %% Solar radiation pressure
    a_SRP = srpAccel(rVec, rSun, const);

    %% Total acceleration
    % a_total = a_grav ...
    %         + a_J2 ...
    %         + a_Sun_grav ...
    %         + a_Moon_grav ...
    %         + a_SRP;
    % 
    % % If you want to include drag too, use this instead:
    % % a_total = a_grav + a_J2 + a_Sun_grav + a_Moon_grav + a_SRP + a_drag;

    %% Total acceleration with switches
    a_total = a_grav;

    if const.useJ2
        a_total = a_total + a_J2;
    end
    
    if const.useSun
        a_total = a_total + a_Sun_grav;
    end
    
    if const.useMoon
        a_total = a_total + a_Moon_grav;
    end
    
    if const.useSRP
        a_total = a_total + a_SRP;
    end
    
    if const.useDrag
        a_total = a_total + a_drag;
    end

    %% Return derivative of state vector
    dydt = [vx; vy; vz; a_total];

end

%% HELPER FUNCTIONS %%
function rSun = sunPositionApprox(t, const)

    thetaS = const.nSun * t;

    rSun = const.AU * [cos(thetaS);
                       sin(thetaS);
                       0];

end

function rMoon = moonPositionApprox(t, const)

    thetaM = const.nMoon * t;

    rMoon = const.rMoonMag * [cos(thetaM);
                              sin(thetaM);
                              0];

end

function a3B = thirdBodyAccel(rSat, rBody, muBody)

    % rSat  = satellite position relative to Earth, km
    % rBody = Sun or Moon position relative to Earth, km
    % muBody = gravitational parameter of Sun or Moon, km^3/s^2

    r_sat_to_body = rBody - rSat;

    direct_term = r_sat_to_body / norm(r_sat_to_body)^3;
    indirect_term = rBody / norm(rBody)^3;

    a3B = muBody * (direct_term - indirect_term);

end

function aSRP = srpAccel(rSat, rSun, const)

    % rSat = satellite position relative to Earth, km
    % rSun = Sun position relative to Earth, km

    % Vector from Sun to satellite
    rSunToSat = rSat - rSun; % km

    % Distance from Sun to satellite
    dSunToSat = norm(rSunToSat); % km

    % Unit vector pointing away from the Sun
    sHat = rSunToSat / dSunToSat;

    % SRP acceleration magnitude in m/s^2
    aSRP_mag_mps2 = const.Psrp * const.Cr * const.A_m * ...
                    (const.AU / dSunToSat)^2;

    % Convert to km/s^2
    aSRP_mag = aSRP_mag_mps2 / 1000;

    % Vector SRP acceleration
    aSRP = aSRP_mag * sHat;

end