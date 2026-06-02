%% MAE 240 Final Project
% GEO / high area-to-mass debris propagation with perturbations
% Model: 2BP + J2 + Sun/Moon third-body gravity + SRP

clc; clear; close all;

%% ============================================================
% 1. Constants
% ============================================================

% Earth parameters
const.mu = 398600.4354;           % Earth gravitational parameter, km^3/s^2
const.R = 6378.1366;              % Earth radius, km
const.J2 = 1.08262668e-3;         % Earth's J2 coefficient
const.omega_E = 7.2921159e-5;     % Earth rotation rate, rad/s

% Third-body gravitational parameters
const.muSun = 1.32712440018e11;   % Sun gravitational parameter, km^3/s^2
const.muMoon = 4902.800066;       % Moon gravitational parameter, km^3/s^2

% Approximate Sun/Moon distances from Earth
const.AU = 149597870.7;           % Astronomical unit, km
const.rMoonMag = 384400;          % Average Earth-Moon distance, km

% Approximate Sun/Moon periods
const.TSun = 365.25 * 24 * 3600;        % apparent Sun period, s
const.TMoon = 27.321661 * 24 * 3600;    % Moon sidereal period, s

const.nSun = 2*pi / const.TSun;         % rad/s
const.nMoon = 2*pi / const.TMoon;       % rad/s

% SRP parameters
const.Psrp = 4.56e-6;             % solar radiation pressure at 1 AU, N/m^2
const.Cr = 1.2;                   % reflectivity coefficient
const.A_m = 0.0123;               % default area-to-mass ratio, m^2/kg

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
const.useDrag = false;  % Keep false for GEO/xGEO cases

%% ============================================================
% 3. Nominal GEO mission parameters
% ============================================================

% GEO-like initial orbit mission parameters
a0 = 42164;          % semi-major axis, km
e0 = 0.001;          % eccentricity
i0 = 0.1;            % inclination, deg
RAAN0 = 0;           % right ascension of ascending node, deg
omega0 = 0;          % argument of perigee, deg
M0 = 0;              % mean anomaly, deg

% Time offset used inside oe_to_cartesian
delta_t = 0;         % seconds

% Convert nominal orbital elements to Cartesian state
X0 = oe_to_cartesian(a0, e0, i0, RAAN0, omega0, M0, const.mu, delta_t);

% Extract initial position and velocity if needed
r0 = X0(1:3);        % km
v0 = X0(4:6);        % km/s

fprintf('Nominal initial radius magnitude: %.3f km\n', norm(r0));
fprintf('Nominal initial velocity magnitude: %.6f km/s\n', norm(v0));

%% ============================================================
% 4. Time span
% ============================================================

days = 364.25;                         % propagation duration, days
t0 = 0;                                % initial time, s
tf = days * 24 * 3600;                 % final time, s

% Output times
Nout = 2000;
tspan = linspace(t0, tf, Nout);
timeDays = tspan(:) / (24*3600);

%% ============================================================
% 5. Numerical integration options
% ============================================================

opts = odeset( ...
    'RelTol', 1e-10, ...
    'AbsTol', 1e-12);

%% ============================================================
% 5A. Two-body validation printout only
% ============================================================

% Purpose:
%   Show that without perturbations, the orbital elements remain constant.
%   This validates the numerical integration and element conversion.
%
%   No plots are made here. Only values are printed.

const2BP = const;

% Turn off all perturbations
const2BP.useJ2 = false;
const2BP.useSun = false;
const2BP.useMoon = false;
const2BP.useSRP = false;
const2BP.useEarthShadow = false;
const2BP.useDrag = false;

[t2BP, X2BP] = ode45(@(t,X) TwoBodyProblem_withPertubations(t, X, const2BP), ...
                     tspan, X0, opts);

[a2BP, e2BP, i2BP, RAAN2BP, omega2BP, M2BP] = ...
    getElementHistory(t2BP, X2BP, const.mu);

energy2BP = zeros(length(t2BP),1);
h2BP = zeros(length(t2BP),1);

for k = 1:length(t2BP)

    r_k = X2BP(k,1:3)';
    v_k = X2BP(k,4:6)';

    energy2BP(k) = 0.5*dot(v_k,v_k) - const.mu/norm(r_k);
    h2BP(k) = norm(cross(r_k,v_k));

end

fprintf('\n============================================================\n');
fprintf('TWO-BODY VALIDATION: ALL PERTURBATIONS OFF\n');
fprintf('============================================================\n');
fprintf('Initial a        = %.12f km\n', a2BP(1));
fprintf('Final a          = %.12f km\n', a2BP(end));
fprintf('Max |a - a0|     = %.6e km\n\n', max(abs(a2BP - a2BP(1))));

fprintf('Initial e        = %.12e\n', e2BP(1));
fprintf('Final e          = %.12e\n', e2BP(end));
fprintf('Max |e - e0|     = %.6e\n\n', max(abs(e2BP - e2BP(1))));

fprintf('Initial i        = %.12e deg\n', i2BP(1));
fprintf('Final i          = %.12e deg\n', i2BP(end));
fprintf('Max |i - i0|     = %.6e deg\n\n', max(abs(i2BP - i2BP(1))));

fprintf('Initial energy   = %.12e km^2/s^2\n', energy2BP(1));
fprintf('Final energy     = %.12e km^2/s^2\n', energy2BP(end));
fprintf('Max |E - E0|     = %.6e km^2/s^2\n\n', max(abs(energy2BP - energy2BP(1))));

fprintf('Initial |h|      = %.12e km^2/s\n', h2BP(1));
fprintf('Final |h|        = %.12e km^2/s\n', h2BP(end));
fprintf('Max |h - h0|     = %.6e km^2/s\n', max(abs(h2BP - h2BP(1))));
fprintf('============================================================\n\n');
%% ============================================================
% 6. Study 1: Fixed GEO orbit, varying area-to-mass ratio
%    Separate figure for each A/m case
% ============================================================

% Fixed GEO mission orbit
a_GEO = 42164;       % km
e_GEO = e0;
i_GEO = i0;
RAAN_GEO = RAAN0;
omega_GEO = omega0;
M_GEO = M0;

% Three representative area-to-mass cases
AmCases = [0.01, 0.05, 10]; % m^2/kg

AmCaseNames = ["Dense satellite, A/m = 0.01 m^2/kg", ...
               "CubeSat / small satellite, A/m = 0.05 m^2/kg", ...
               "HAMR debris, A/m = 10 m^2/kg"];

% Save original A/m
A_m_original = const.A_m;

for j = 1:length(AmCases)

    const.A_m = AmCases(j);

    % Same GEO initial condition for every A/m case
    X0_case = oe_to_cartesian(a_GEO, e_GEO, i_GEO, RAAN_GEO, ...
                              omega_GEO, M_GEO, const.mu, delta_t);

    [tOut_case, XOut_case] = ode45(@(t,X) TwoBodyProblem_withPertubations(t, X, const), ...
                                   tspan, X0_case, opts);

    [aHist, eHist, iHist, RAANHist, omegaHist, MHist] = ...
        getElementHistory(tOut_case, XOut_case, const.mu);

    timeDays_case = tOut_case / (24*3600);

    plotSixOE_singleCase(timeDays_case, aHist, eHist, iHist, RAANHist, omegaHist, MHist, ...
                         "GEO Orbit: " + AmCaseNames(j));

end

% Restore original A/m
const.A_m = A_m_original;


%% ============================================================
% 7. Study 2: Fixed area-to-mass ratio, varying altitude
%    Separate figure for each altitude case
% ============================================================

% Fixed satellite area-to-mass ratio
const.A_m = 10; % m^2/kg, large debris

% Three high-altitude orbit cases above Earth's surface
% Drag is negligible for all three.
altCases = [35786, 36086, 100000]; % km

orbitNames = ["GEO, h = 35786 km", ...
              "Graveyard, h = GEO + 300 km", ...
              "High Earth Orbit, h = 100000 km"];

% Convert altitude above Earth to semi-major axis
aCases = const.R + altCases;

for j = 1:length(aCases)

    a0_case = aCases(j);

    % Same eccentricity/orientation parameters, different distance from Earth
    X0_case = oe_to_cartesian(a0_case, e0, i0, RAAN0, ...
                              omega0, M0, const.mu, delta_t);

    [tOut_case, XOut_case] = ode45(@(t,X) TwoBodyProblem_withPertubations(t, X, const), ...
                                   tspan, X0_case, opts);

    [aHist, eHist, iHist, RAANHist, omegaHist, MHist] = ...
        getElementHistory(tOut_case, XOut_case, const.mu);

    timeDays_case = tOut_case / (24*3600);

    plotSixOE_singleCase(timeDays_case, aHist, eHist, iHist, RAANHist, omegaHist, MHist, ...
                         "Fixed A/m = 10 m^2/kg: " + orbitNames(j));

end

% Restore original A/m
const.A_m = A_m_original;

%% ============================================================
% 8. Optional: nominal GEO trajectory and Earth shadow fraction
% ============================================================

const.A_m = A_m_original;

X0_nominal = oe_to_cartesian(a0, e0, i0, RAAN0, omega0, M0, const.mu, delta_t);

[tNom, XNom] = ode45(@(t,X) TwoBodyProblem_withPertubations(t, X, const), ...
                     tspan, X0_nominal, opts);

rNom = XNom(:,1:3);

figure('Name', 'Nominal GEO Trajectory', 'Color', 'w');
plot3(rNom(:,1), rNom(:,2), rNom(:,3), 'LineWidth', 1.5);
hold on;
plot3(0, 0, 0, 'ko', 'MarkerFaceColor', 'k');
grid on; axis equal;
xlabel('x [km]');
ylabel('y [km]');
zlabel('z [km]');
title('Nominal GEO Trajectory in Earth-Centered Inertial Frame');

shadowFraction = computeShadowFraction(tNom, XNom, const);
fprintf('Fraction of nominal trajectory in Earth shadow: %.4f%%\n', shadowFraction*100);

%% ============================================================
% 9. Optional: nominal two-body-style validation quantities
% ============================================================
% Note: With perturbations ON, energy and angular momentum are not expected
% to remain exactly constant. For pure 2BP validation, turn all perturbation
% switches off and rerun this section.

rHist = XNom(:,1:3);
vHist = XNom(:,4:6);

energyHist = zeros(length(tNom),1);
hHist = zeros(length(tNom),1);

for k = 1:length(tNom)

    r_k = rHist(k,:)';
    v_k = vHist(k,:)';

    energyHist(k) = 0.5*dot(v_k,v_k) - const.mu/norm(r_k);
    hHist(k) = norm(cross(r_k,v_k));

end

figure('Name', 'Specific Energy Change', 'Color', 'w');
plot(tNom/(24*3600), energyHist - energyHist(1), 'LineWidth', 1.5);
grid on;
xlabel('Time [days]');
ylabel('\Delta specific energy [km^2/s^2]');
title('Specific Energy Change');

figure('Name', 'Angular Momentum Change', 'Color', 'w');
plot(tNom/(24*3600), hHist - hHist(1), 'LineWidth', 1.5);
grid on;
xlabel('Time [days]');
ylabel('\Delta |h| [km^2/s]');
title('Angular Momentum Magnitude Change');

%% ============================================================
% 10. Fourteen-day CubeSat drift from target orbit
% ============================================================

% Purpose:
%   This section does NOT implement control.
%   It samples the uncontrolled perturbed CubeSat orbit every 14 days and
%   compares the osculating elements to the target GEO orbit.
%
%   These 14-day errors show what future correction/control would need to
%   reduce.

% CubeSat / small satellite case already used in Study 1
constCubeSat = const;
constCubeSat.A_m = 0.05;  % m^2/kg

% Target GEO orbit
aTarget = a0;
eTarget = e0;
iTarget = i0;

% Initial condition is the target GEO orbit
X0_cube = oe_to_cartesian(aTarget, eTarget, iTarget, ...
                          RAAN0, omega0, M0, const.mu, delta_t);

% Propagate CubeSat case with perturbations ON
[tCube, XCube] = ode45(@(t,X) TwoBodyProblem_withPertubations(t, X, constCubeSat), ...
                       tspan, X0_cube, opts);

[aCube, eCube, iCube, RAANCube, omegaCube, MCube] = ...
    getElementHistory(tCube, XCube, const.mu);

timeDaysCube = tCube / (24*3600);

% Sample every 14 days
sampleIntervalDays = 14;
sampleDays = 0:sampleIntervalDays:days;

aSample = interp1(timeDaysCube, aCube, sampleDays, 'linear');
eSample = interp1(timeDaysCube, eCube, sampleDays, 'linear');
iSample = interp1(timeDaysCube, iCube, sampleDays, 'linear');

% Difference between perturbed osculating elements and target orbit
deltaA_sample = aSample - aTarget;
deltaE_sample = eSample - eTarget;
deltaI_sample = iSample - iTarget;

% Print summary values
fprintf('\n============================================================\n');
fprintf('14-DAY CUBESAT DRIFT FROM TARGET ORBIT\n');
fprintf('============================================================\n');
fprintf('CubeSat A/m                    = %.4f m^2/kg\n', constCubeSat.A_m);
fprintf('Max sampled |a - aTarget|      = %.6f km\n', max(abs(deltaA_sample)));
fprintf('Max sampled |e - eTarget|      = %.6e\n', max(abs(deltaE_sample)));
fprintf('Max sampled |i - iTarget|      = %.6e deg\n', max(abs(deltaI_sample)));
fprintf('Final sampled a - aTarget      = %.6f km\n', deltaA_sample(end));
fprintf('Final sampled e - eTarget      = %.6e\n', deltaE_sample(end));
fprintf('Final sampled i - iTarget      = %.6e deg\n', deltaI_sample(end));
fprintf('============================================================\n\n');

% One final graph showing the 14-day drift in a, e, and i
figure('Name', '14-Day CubeSat Drift from Target Orbit', 'Color', 'w');
tiledlayout(3,1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
plot(sampleDays, deltaA_sample, 'o-', 'LineWidth', 1.5);
grid on;
ylabel('a - a_{target} [km]');
title('14-Day CubeSat Drift from Target Orbit');

nexttile;
plot(sampleDays, deltaE_sample, 'o-', 'LineWidth', 1.5);
grid on;
ylabel('e - e_{target}');

nexttile;
plot(sampleDays, deltaI_sample, 'o-', 'LineWidth', 1.5);
grid on;
xlabel('Time [days]');
ylabel('i - i_{target} [deg]');

%% ============================================================
% Local functions
% ============================================================

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

function plotSixOE_singleCase(timeDays, aHist, eHist, iHist, RAANHist, omegaHist, MHist, ...
                              mainTitleText)

    figure('Name', mainTitleText, 'Color', 'w');

    tiledlayout(3,2, 'TileSpacing', 'compact', 'Padding', 'compact');

    % Semi-major axis
    nexttile;
    plot(timeDays, aHist, 'LineWidth', 1.4);
    grid on;
    xlabel('Time [days]');
    ylabel('a [km]');
    title('Semi-major Axis');

    % Eccentricity
    nexttile;
    plot(timeDays, eHist, 'LineWidth', 1.4);
    grid on;
    xlabel('Time [days]');
    ylabel('e');
    title('Eccentricity');

    % Inclination
    nexttile;
    plot(timeDays, iHist, 'LineWidth', 1.4);
    grid on;
    xlabel('Time [days]');
    ylabel('i [deg]');
    title('Inclination');

    % RAAN
    nexttile;
    plot(timeDays, RAANHist, 'LineWidth', 1.4);
    grid on;
    xlabel('Time [days]');
    ylabel('\Omega [deg]');
    title('RAAN');

    % Argument of perigee
    nexttile;
    plot(timeDays, omegaHist, 'LineWidth', 1.4);
    grid on;
    xlabel('Time [days]');
    ylabel('\omega [deg]');
    title('Argument of Perigee');

    % Mean anomaly
    nexttile;
    plot(timeDays, MHist, 'LineWidth', 1.4);
    grid on;
    xlabel('Time [days]');
    ylabel('M [deg]');
    title('Mean Anomaly');

    sgtitle(mainTitleText);

end

function shadowFraction = computeShadowFraction(tOut, XOut, const)

    inShadowArray = false(length(tOut), 1);

    for k = 1:length(tOut)

        t = tOut(k);
        rSat = XOut(k,1:3)';

        rSun = sunPositionApprox_main(t, const);

        inShadowArray(k) = earthShadowCylindrical_main(rSat, rSun, const);

    end

    shadowFraction = sum(inShadowArray) / length(inShadowArray);

end

function inShadow = earthShadowCylindrical_main(rSat, rSun, const)

    % Unit vector from Earth to Sun
    sHat_EarthToSun = rSun / norm(rSun);

    % Projection of satellite position onto Sun direction
    projection = dot(rSat, sHat_EarthToSun);

    % Perpendicular distance from Sun-Earth line
    rPerp = norm(rSat - projection * sHat_EarthToSun);

    % Cylindrical shadow approximation
    inShadow = (projection < 0) && (rPerp < const.R);

end

function rSun = sunPositionApprox_main(t, const)

    thetaS = const.nSun * t;

    rSun = const.AU * [cos(thetaS);
                       sin(thetaS);
                       0];

end