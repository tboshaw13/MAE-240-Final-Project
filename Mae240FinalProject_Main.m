%% MAE 240 Final Project
% GEO / high area-to-mass debris propagation with perturbations
% Model: 2BP + J2 + Sun/Moon third-body gravity + SRP

clc; clear; close all;

%% ============================================================
% 1. Constants
% ============================================================

% Earth parameters
const.mu = 398600.4354;          % Earth gravitational parameter, km^3/s^2
const.R = 6378.1366;             % Earth radius, km
const.J2 = 1.08262668e-3;        % Earth's J2 coefficient
const.omega_E = 7.2921159e-5;    % Earth rotation rate, rad/s

% Third-body gravitational parameters
const.muSun = 1.32712440018e11;  % Sun gravitational parameter, km^3/s^2
const.muMoon = 4902.800066;      % Moon gravitational parameter, km^3/s^2

% Approximate Sun/Moon distances from Earth
const.AU = 149597870.7;          % Astronomical unit, km
const.rMoonMag = 384400;         % Average Earth-Moon distance, km

% Approximate Sun/Moon periods
const.TSun = 365.25 * 24 * 3600;       % apparent Sun period, s
const.TMoon = 27.321661 * 24 * 3600;   % Moon sidereal period, s

const.nSun = 2*pi / const.TSun;        % rad/s
const.nMoon = 2*pi / const.TMoon;      % rad/s

% SRP parameters
const.Psrp = 4.56e-6;            % solar radiation pressure at 1 AU, N/m^2
const.Cr = 1.2;                  % reflectivity coefficient
const.A_m = 0.0123;              % area-to-mass ratio, m^2/kg

% Drag parameters, optional
const.Cd = 1.28;

%% ============================================================
% 2. Model switches
% ============================================================

const.useJ2 = true;
const.useSun = true;
const.useMoon = true;
const.useSRP = true;
const.useEarthShadow = true;
const.useDrag = false;  % Keep false for GEO cases

%% ============================================================
% 3. Initial mission orbit using orbital elements
% ============================================================

% GEO-like initial orbit
a0 = 42164;          % semi-major axis, km
e0 = 0.001;          % eccentricity
i0 = 0.1;            % inclination, deg
RAAN0 = 0;           % right ascension of ascending node, deg
omega0 = 0;          % argument of perigee, deg
M0 = 0;              % mean anomaly, deg

% Time offset used inside oe_to_cartesian
delta_t = 0;         % seconds

% Convert orbital elements to Cartesian state
X0 = oe_to_cartesian(a0, e0, i0, RAAN0, omega0, M0, const.mu, delta_t);

% Extract initial position and velocity if needed
r0 = X0(1:3);        % km
v0 = X0(4:6);        % km/s

%% ============================================================
% 4. Time span
% ============================================================

days = 364.25;                         % propagation duration
t0 = 0;
tf = days * 24 * 3600;             % seconds

% Output times. This helps give smooth plots.
Nout = 2000;
tspan = linspace(t0, tf, Nout);

%% ============================================================
% 5. Numerical integration with ode45
% ============================================================

opts = odeset( ...
    'RelTol', 1e-10, ...
    'AbsTol', 1e-12);

[tOut, XOut] = ode45(@(t,X) TwoBodyProblem_withPertubations(t, X, const), ...
                     tspan, X0, opts);

%% ============================================================
% 6. Extract r and v histories
% ============================================================

rHist = XOut(:,1:3);   % km
vHist = XOut(:,4:6);   % km/s

%% ============================================================
% 7. Convert r,v history back to orbital elements
% ============================================================

aHist = zeros(length(tOut), 1);
eHist = zeros(length(tOut), 1);
iHist = zeros(length(tOut), 1);
RAANHist = zeros(length(tOut), 1);
omegaHist = zeros(length(tOut), 1);
MHist = zeros(length(tOut), 1);

% Important:
% For each time step, use delta_t = 0 here.
% The integrated Cartesian state already corresponds to that time.
% We are only converting that instantaneous state into osculating elements.
delta_t_oe = 0;

for k = 1:length(tOut)

    X_k = XOut(k,:)';

    [aHist(k), eHist(k), iHist(k), RAANHist(k), omegaHist(k), MHist(k)] = ...
        Cartesian_to_oe(X_k, const.mu, delta_t_oe);

end

%% ============================================================
% 8. Basic plots
% ============================================================

timeDays = tOut / (24*3600);

figure;
plot(timeDays, aHist, 'LineWidth', 1.5);
grid on;
xlabel('Time [days]');
ylabel('Semi-major axis a [km]');
title('Osculating Semi-major Axis');

figure;
plot(timeDays, eHist, 'LineWidth', 1.5);
grid on;
xlabel('Time [days]');
ylabel('Eccentricity e');
title('Osculating Eccentricity');

figure;
plot(timeDays, iHist, 'LineWidth', 1.5);
grid on;
xlabel('Time [days]');
ylabel('Inclination i [deg]');
title('Osculating Inclination');

figure;
plot(timeDays, RAANHist, 'LineWidth', 1.5);
grid on;
xlabel('Time [days]');
ylabel('RAAN \Omega [deg]');
title('Osculating RAAN');

figure;
plot(timeDays, omegaHist, 'LineWidth', 1.5);
grid on;
xlabel('Time [days]');
ylabel('Argument of Perigee \omega [deg]');
title('Osculating Argument of Perigee');

figure;
plot(timeDays, MHist, 'LineWidth', 1.5);
grid on;
xlabel('Time [days]');
ylabel('Mean Anomaly M [deg]');
title('Osculating Mean Anomaly');

%% ============================================================
% 9. Two-body validation quantities
% ============================================================

energyHist = zeros(length(tOut),1);
hHist = zeros(length(tOut),1);

for k = 1:length(tOut)

    r_k = rHist(k,:)';
    v_k = vHist(k,:)';

    energyHist(k) = 0.5*dot(v_k,v_k) - const.mu/norm(r_k);
    hHist(k) = norm(cross(r_k,v_k));

end

figure;
plot(timeDays, energyHist - energyHist(1), 'LineWidth', 1.5);
grid on;
xlabel('Time [days]');
ylabel('\Delta specific energy [km^2/s^2]');
title('Specific Energy Change');

figure;
plot(timeDays, hHist - hHist(1), 'LineWidth', 1.5);
grid on;
xlabel('Time [days]');
ylabel('\Delta |h| [km^2/s]');
title('Angular Momentum Magnitude Change');

%% ============================================================
% 10. Area-to-mass ratio sweep
% ============================================================

% Save original value so we can restore it later
A_m_original = const.A_m;

% Area-to-mass ratio cases, m^2/kg
AmCases = [0.005, 0.0123, 0.05, 0.1, 0.5, 1.0];

maxE_Am = zeros(length(AmCases), 1);
deltaA_Am = zeros(length(AmCases), 1);
deltaI_Am = zeros(length(AmCases), 1);

% Use fewer output points for sweeps so the code runs faster
NoutSweep = 800;
tspanSweep = linspace(t0, tf, NoutSweep);

for j = 1:length(AmCases)

    const.A_m = AmCases(j);

    % Same initial orbit for every A/m case
    X0_sweep = oe_to_cartesian(a0, e0, i0, RAAN0, omega0, M0, const.mu, delta_t);

    [tSweep, XSweep] = ode45(@(t,X) TwoBodyProblem_withPertubations(t, X, const), ...
                             tspanSweep, X0_sweep, opts);

    [aSweep, eSweep, iSweep, ~, ~, ~] = ...
        getElementHistory(tSweep, XSweep, const.mu);

    maxE_Am(j) = max(eSweep);
    deltaA_Am(j) = aSweep(end) - aSweep(1);
    deltaI_Am(j) = iSweep(end) - iSweep(1);

end

% Restore original A/m
const.A_m = A_m_original;

figure;
plot(AmCases, maxE_Am, 'o-', 'LineWidth', 1.5);
grid on;
xlabel('Area-to-mass ratio A/m [m^2/kg]');
ylabel('Maximum eccentricity');
title('Effect of Area-to-Mass Ratio on GEO Eccentricity');

figure;
plot(AmCases, deltaA_Am, 'o-', 'LineWidth', 1.5);
grid on;
xlabel('Area-to-mass ratio A/m [m^2/kg]');
ylabel('\Delta a over propagation [km]');
title('Effect of Area-to-Mass Ratio on Semi-major Axis Drift');

figure;
plot(AmCases, deltaI_Am, 'o-', 'LineWidth', 1.5);
grid on;
xlabel('Area-to-mass ratio A/m [m^2/kg]');
ylabel('\Delta i over propagation [deg]');
title('Effect of Area-to-Mass Ratio on Inclination Drift');

%% ============================================================
% 11. Distance-from-Earth / altitude sweep
% ============================================================

% Initial altitude cases above Earth surface, km
% GEO altitude is about 35786 km.
altCases = [35786, 37000, 40000, 45000, 50000];

aCases = const.R + altCases;

maxE_alt = zeros(length(aCases), 1);
deltaA_alt = zeros(length(aCases), 1);
deltaI_alt = zeros(length(aCases), 1);

% Keep original A/m for this sweep
const.A_m = A_m_original;

for ia = 1:length(aCases)

    a0_case = aCases(ia);

    X0_sweep = oe_to_cartesian(a0_case, e0, i0, RAAN0, omega0, M0, ...
                               const.mu, delta_t);

    [tSweep, XSweep] = ode45(@(t,X) TwoBodyProblem_withPertubations(t, X, const), ...
                             tspanSweep, X0_sweep, opts);

    [aSweep, eSweep, iSweep, ~, ~, ~] = ...
        getElementHistory(tSweep, XSweep, const.mu);

    maxE_alt(ia) = max(eSweep);
    deltaA_alt(ia) = aSweep(end) - aSweep(1);
    deltaI_alt(ia) = iSweep(end) - iSweep(1);

end

figure;
plot(altCases, maxE_alt, 'o-', 'LineWidth', 1.5);
grid on;
xlabel('Initial altitude above Earth [km]');
ylabel('Maximum eccentricity');
title('Effect of Distance from Earth on Eccentricity');

figure;
plot(altCases, deltaA_alt, 'o-', 'LineWidth', 1.5);
grid on;
xlabel('Initial altitude above Earth [km]');
ylabel('\Delta a over propagation [km]');
title('Effect of Distance from Earth on Semi-major Axis Drift');

figure;
plot(altCases, deltaI_alt, 'o-', 'LineWidth', 1.5);
grid on;
xlabel('Initial altitude above Earth [km]');
ylabel('\Delta i over propagation [deg]');
title('Effect of Distance from Earth on Inclination Drift');

%% ============================================================
% 12. Heatmap: altitude and area-to-mass ratio sweep
% ============================================================

% Matrix rows = altitude cases
% Matrix columns = A/m cases
maxEGrid = zeros(length(altCases), length(AmCases));
deltaAGrid = zeros(length(altCases), length(AmCases));

for ia = 1:length(aCases)

    for j = 1:length(AmCases)

        a0_case = aCases(ia);
        const.A_m = AmCases(j);

        X0_sweep = oe_to_cartesian(a0_case, e0, i0, RAAN0, omega0, M0, ...
                                   const.mu, delta_t);

        [tSweep, XSweep] = ode45(@(t,X) TwoBodyProblem_withPertubations(t, X, const), ...
                                 tspanSweep, X0_sweep, opts);

        [aSweep, eSweep, ~, ~, ~, ~] = ...
            getElementHistory(tSweep, XSweep, const.mu);

        maxEGrid(ia,j) = max(eSweep);
        deltaAGrid(ia,j) = aSweep(end) - aSweep(1);

    end

end

% Restore original value
const.A_m = A_m_original;

figure;
imagesc(AmCases, altCases, maxEGrid);
set(gca, 'YDir', 'normal');
colorbar;
xlabel('Area-to-mass ratio A/m [m^2/kg]');
ylabel('Initial altitude above Earth [km]');
title('Maximum Eccentricity vs Altitude and Area-to-Mass Ratio');

figure;
imagesc(AmCases, altCases, deltaAGrid);
set(gca, 'YDir', 'normal');
colorbar;
xlabel('Area-to-mass ratio A/m [m^2/kg]');
ylabel('Initial altitude above Earth [km]');
title('\Delta a vs Altitude and Area-to-Mass Ratio');

shadowFraction = computeShadowFraction(tOut, XOut, const);
fprintf('Fraction of trajectory in Earth shadow: %.4f%%\n', shadowFraction*100);



%% Functions
function [aHist, eHist, iHist, RAANHist, omegaHist, MHist] = ...
    getElementHistory(tOut, XOut, mu)

    aHist = zeros(length(tOut), 1);
    eHist = zeros(length(tOut), 1);
    iHist = zeros(length(tOut), 1);
    RAANHist = zeros(length(tOut), 1);
    omegaHist = zeros(length(tOut), 1);
    MHist = zeros(length(tOut), 1);

    delta_t_oe = 0;

    for k = 1:length(tOut)

        X_k = XOut(k,:)';

        [aHist(k), eHist(k), iHist(k), RAANHist(k), omegaHist(k), MHist(k)] = ...
            Cartesian_to_oe(X_k, mu, delta_t_oe);

    end

end

function shadowFraction = computeShadowFraction(tOut, XOut, const)

    inShadowArray = false(length(tOut), 1);

    for k = 1:length(tOut)

        t = tOut(k);
        rSat = XOut(k,1:3)';

        rSun = sunPositionApprox(t, const);

        inShadowArray(k) = earthShadowCylindrical_main(rSat, rSun, const);

    end

    shadowFraction = sum(inShadowArray) / length(inShadowArray);

end


function inShadow = earthShadowCylindrical_main(rSat, rSun, const)

    sHat_EarthToSun = rSun / norm(rSun);

    projection = dot(rSat, sHat_EarthToSun);

    rPerp = norm(rSat - projection * sHat_EarthToSun);

    inShadow = (projection < 0) && (rPerp < const.R);

end

function rSun = sunPositionApprox(t, const)

    thetaS = const.nSun * t;

    rSun = const.AU * [cos(thetaS);
                       sin(thetaS);
                       0];

end