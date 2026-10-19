%% DR2, Coverage Planning: kinematic lawnmower model + Monte Carlo
% Advances the DR1 coverage trade (Trade Studies/Coverage Path Planning Trade)
% by adding: the camera's operating resolution/FOV, trapezoidal accel/decel on
% every leg, stop-and-settle turns, path-following + navigation error,
% attitude-induced footprint shift, discrete image sampling with motion blur,
% and a mission time / energy budget.
%
% Geometric coverage (area seen by any frame) and marker recognition (marker
% fully in a frame with enough pixels and low enough blur) are kept separate.
%
% Values tagged TBD are placeholders to be replaced with team/datasheet values.

%% Housekeeping
clc, clear, close all;
rng(1);

scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir), scriptDir = pwd; end
figDir = fullfile(scriptDir, 'figures');
if ~exist(figDir, 'dir'), mkdir(figDir); end

%% Monte Carlo settings
nRuns = 100;
nRunsSweep = 30;

%% Search area and flight envelope (GNC requirements draft)
p.areaWidth = 5;              % [m]
p.areaLength = 5;             % [m]
p.altitude = 2.0;             % [m] AGL search altitude (h_max = 2.0 m)
p.vCruise = 0.40;             % [m/s] nominal cruise speed
p.vMax = 0.50;                % [m/s] maximum flight speed
p.tiltLimit = 10;             % [deg] nominal roll/pitch limit
p.g = 9.81;                   % [m/s^2]

%% Path-following assumptions (TBD: confirm with controls / localization)
p.accel = 0.50;               % [m/s^2] commanded accel/decel on each leg (TBD)
p.cornerHold = 1.0;           % [s] stop-and-settle at each waypoint (TBD)
p.overlap = 0.20;             % [-] minimum overlap between adjacent tracks
p.posErrorStd = 0.05;         % [m] 1-sigma horizontal error, tracking + nav, per axis (TBD)
p.altitudeErrorStd = 0.05;    % [m] 1-sigma altitude error (TBD)
p.attitudeJitterStd = 1.0;    % [deg] 1-sigma roll/pitch jitter about commanded tilt (TBD)
p.posErrorTau = 2.0;          % [s] correlation time of position/altitude error
p.attitudeTau = 0.5;          % [s] correlation time of attitude jitter

%% Camera (TBD: replace with selected camera datasheet values)
p.imageResolution = [1920, 1080]; % [px] operating resolution (DR1: >= 1080p)
p.hFOV = 60;                  % [deg] horizontal FOV, long axis mounted cross-track
p.frameRate = 10;             % [Hz] frames processed for detection on the Pi 5 (TBD)
p.exposure = 0.005;           % [s] exposure time (TBD)
p.vFOV = 2*atand((p.imageResolution(2)/p.imageResolution(1))*tand(p.hFOV/2));

%% Marker recognition criteria
p.markerSize = 0.08;          % [m] ArUco side length
p.minimumMarkerPixels = 60;   % [px] per side for reliable decode (GSD trade)
p.maxBlurPixels = 2.0;        % [px] maximum motion blur for decode (TBD)
p.minRecognitionFrames = 3;   % [-] good frames needed to confirm the marker (TBD)

%% Analysis grids
p.gridSpacing = 0.05;         % [m] geometric coverage grid
p.markerGridSpacing = 0.10;   % [m] candidate marker locations
p.dt = 0.01;                  % [s] trajectory time step

%% Mission timeline (TBD: confirm with ConOps)
m.launchPoint = [-2, 0];      % [m] launch site relative to search-area corner (TBD)
m.initHover = 5;              % [s] initialization hover before transit (TBD)
m.takeoffTime = 3.0;          % [s] climb to search altitude (T_s,takeoff <= 3.0 s)
m.markerHover = 5;            % [s] hover to confirm marker (TBD)
m.commsHover = 5;             % [s] hover while coordinates are transmitted (TBD)
m.descentRate = 0.3;          % [m/s] landing descent rate (TBD)

%% Energy inputs (TBD: from power/propulsion)
e.mass = 1.5;                 % [kg] all-up mass (TBD)
e.hoverPower = 180;           % [W] electrical hover power (TBD)
e.climbEfficiency = 0.5;      % [-] electrical-to-climb efficiency (TBD)
e.batteryEnergy = 74;         % [Wh] e.g. 4S 5000 mAh (TBD)
e.usableFraction = 0.80;      % [-] usable depth of discharge (TBD)
e.reserveFraction = 0.30;     % [-] energy reserve kept at landing (TBD)

%% Requirements (TBD: replace with approved values)
req.coverage = 100;           % [%] geometric coverage of the search area
req.recognition = 99;         % [%] marker locations recognized (TBD)
req.maxGSD = p.markerSize/p.minimumMarkerPixels; % [m/px] derived from pixel requirement
req.searchTime = 120;         % [s] (TBD)
req.missionTime = 300;        % [s] (TBD)
req.energy = e.batteryEnergy*e.usableFraction*(1 - e.reserveFraction); % [Wh] allowed use

%% Footprint and error budget
fp = cameraFootprint(p.altitude, p);

% 3-sigma margins the inset path keeps between the footprint edge and the area
crossTrackMargin = 3*p.posErrorStd + p.altitude*tand(3*p.attitudeJitterStd);
alongTrackMargin = crossTrackMargin + p.altitude*(p.accel/p.g);

fprintf('Footprint at %.1f m: %.2f m x %.2f m, GSD %.2f mm/px, marker %.1f px\n', ...
    p.altitude, fp.width, fp.length, 1000*fp.gsd, p.markerSize/fp.gsd);
fprintf('Path margins: cross-track %.2f m, along-track %.2f m\n', crossTrackMargin, alongTrackMargin);

%% Build commanded paths: DR1 layout vs DR2 inset layout
wpDR1 = generateLawnmowerDR1(p, fp);
wpDR2 = generateInsetLawnmower(p, fp, crossTrackMargin, alongTrackMargin);

% Search starts where the drone reaches the area: the border point closest to launch
entryPoint = min(max(m.launchPoint, [0 0]), [p.areaWidth p.areaLength]);
wpDR2 = [entryPoint; wpDR2];

% Search ends where the drone leaves the area: the border point closest to the last track end
lastPoint = wpDR2(end,:);
edges = [0, p.areaWidth, 0, p.areaLength];
[~, k] = min(abs(lastPoint([1 1 2 2]) - edges));
exitPoint = lastPoint;
exitPoint(ceil(k/2)) = edges(k);
wpDR2 = [wpDR2; exitPoint];

trajDR1 = buildTrajectory(wpDR1, p.vCruise, p.accel, p.cornerHold, p.dt);
trajDR2 = buildTrajectory(wpDR2, p.vCruise, p.accel, p.cornerHold, p.dt);

pathLengthDR1 = sum(vecnorm(diff(wpDR1), 2, 2));
pathLengthDR2 = sum(vecnorm(diff(wpDR2), 2, 2));
dr1ReportedTime = pathLengthDR1/0.5; % DR1 model: constant 0.5 m/s, instant turns

fprintf('\nDR1 layout: %.1f m, DR1-reported %.1f s, kinematic %.1f s\n', ...
    pathLengthDR1, dr1ReportedTime, trajDR1.t(end));
fprintf('DR2 layout: %.1f m, kinematic %.1f s (%d tracks)\n', ...
    pathLengthDR2, trajDR2.t(end), (size(wpDR2,1) - 2)/2);

%% Monte Carlo on the DR2 path
[X, Y] = meshgrid(0:p.gridSpacing:p.areaWidth, 0:p.gridSpacing:p.areaLength);
half = p.markerSize/2;
[MX, MY] = meshgrid(half:p.markerGridSpacing:p.areaWidth-half, half:p.markerGridSpacing:p.areaLength-half);
markerXY = [MX(:), MY(:)];

mc = runMonteCarlo(trajDR2, p, nRuns, X, Y, markerXY);
mcDR1 = runMonteCarlo(trajDR1, p, nRuns, X, Y, markerXY);

%% Mission time and energy budget (worst case: marker found at end of search)
budget = missionBudget(trajDR2, wpDR2, p, m, e);

%% Results vs requirements
fprintf('\n=== DR2 coverage results, %d-run Monte Carlo (DR2 inset path) ===\n', nRuns);
fprintf('%-34s %-14s %-14s %-8s\n', 'Metric', 'Requirement', 'Predicted', 'Status');
printRow('Geometric coverage, worst run [%]', sprintf('>= %.0f', req.coverage), min(mc.coverage), min(mc.coverage) >= req.coverage);
printRow('Marker recognition, worst run [%]', sprintf('>= %.0f', req.recognition), min(mc.recognition), min(mc.recognition) >= req.recognition);
printRow('GSD, worst frame [mm/px]', sprintf('<= %.2f', 1000*req.maxGSD), 1000*max(mc.worstGSD), max(mc.worstGSD) <= req.maxGSD);
printRow('Marker pixels, worst frame [px]', sprintf('>= %d', p.minimumMarkerPixels), min(mc.minPixels), min(mc.minPixels) >= p.minimumMarkerPixels);
printRow('Motion blur, worst frame [px]', sprintf('<= %.1f', p.maxBlurPixels), max(mc.maxBlur), max(mc.maxBlur) <= p.maxBlurPixels);
printRow('Path length [m]', '-', pathLengthDR2, true);
printRow('Search time [s]', sprintf('<= %.0f', req.searchTime), trajDR2.t(end), trajDR2.t(end) <= req.searchTime);
printRow('Worst time to recognition [s]', '-', max(mc.worstRecognitionTime), true);
printRow('Total mission time [s]', sprintf('<= %.0f', req.missionTime), budget.totalTime, budget.totalTime <= req.missionTime);
printRow('Mission energy [Wh]', sprintf('<= %.1f', req.energy), budget.totalEnergy, budget.totalEnergy <= req.energy);

fprintf('\nDR1 layout for comparison: worst coverage %.1f %%, worst recognition %.1f %%\n', ...
    min(mcDR1.coverage), min(mcDR1.recognition));

fprintf('\nMission budget by phase:\n');
for i = 1:numel(budget.phases)
    fprintf('  %-16s %6.1f s  %5.2f Wh\n', budget.phases{i}, budget.time(i), budget.energy(i));
end

%% Sensitivity 1: position error (path fixed, designed for nominal error)
errorValues = [0, 0.05, 0.10, 0.15, 0.20, 0.30];
sweepCoverage = zeros(size(errorValues));
sweepRecognition = zeros(size(errorValues));
for i = 1:numel(errorValues)
    ps = p;
    ps.posErrorStd = errorValues(i);
    mcs = runMonteCarlo(trajDR2, ps, nRunsSweep, X, Y, markerXY);
    sweepCoverage(i) = min(mcs.coverage);
    sweepRecognition(i) = min(mcs.recognition);
end

%% Sensitivity 2: commanded acceleration (search time, analytic)
accelValues = [0.25, 0.5, 1.0, p.g*tand(p.tiltLimit)];
accelSearchTime = zeros(size(accelValues));
for i = 1:numel(accelValues)
    trajA = buildTrajectory(wpDR2, p.vCruise, accelValues(i), p.cornerHold, p.dt);
    accelSearchTime(i) = trajA.t(end);
end

%% Sensitivity 3: marker pixels vs altitude for candidate resolutions
altitudes = linspace(0.5, 2.5, 81);
resolutions = [1280, 720; 1920, 1080; 2560, 1440; 3840, 2160];
markerPixelsVsAltitude = zeros(size(resolutions,1), numel(altitudes));
for i = 1:size(resolutions,1)
    markerPixelsVsAltitude(i,:) = p.markerSize*resolutions(i,1)./(2*altitudes*tand(p.hFOV/2));
end

%% Plots
% Paths and coverage
newFigure([100 100 1400 450]);
subplot(1,3,1);
plotPath(wpDR1, p, sprintf('DR1 layout: %.1f m, %.1f s', pathLengthDR1, trajDR1.t(end)));
subplot(1,3,2);
plotPath(wpDR2, p, sprintf('DR2 inset layout: %.1f m, %.1f s', pathLengthDR2, trajDR2.t(end)));
subplot(1,3,3);
imagesc(X(1,:), Y(:,1), mc.exampleCountGrid); axis xy equal tight; hold on;
plot(mc.exampleCenters(:,1), mc.exampleCenters(:,2), 'w-', 'LineWidth', 0.8);
cb = colorbar; cb.Label.String = 'Frames seeing cell';
xlabel('X Position (m)'); ylabel('Y Position (m)');
title('DR2 coverage map, one Monte Carlo run');
exportgraphics(gcf, fullfile(figDir, 'coverage_paths.png'), 'Resolution', 150);

% Speed profile
newFigure([100 100 900 350]);
plot(trajDR2.t, vecnorm(trajDR2.vel, 2, 2), 'LineWidth', 1.5); hold on;
plot([0 dr1ReportedTime], [0.5 0.5], '--', 'LineWidth', 1.5);
yline(p.vMax, 'k:', 'V_{max}');
xlabel('Time (s)'); ylabel('Speed (m/s)'); ylim([0 0.6]); grid on;
legend('DR2 kinematic model', 'DR1 assumption (constant 0.5 m/s)', 'Location', 'southoutside', 'Orientation', 'horizontal');
title('Search speed profile: accel/decel and stop-and-settle at each waypoint');
exportgraphics(gcf, fullfile(figDir, 'speed_profile.png'), 'Resolution', 150);

% Sensitivities
newFigure([100 100 1400 400]);
subplot(1,3,1);
plot(100*errorValues, sweepCoverage, '-o', 100*errorValues, sweepRecognition, '-s', 'LineWidth', 1.5);
yline(req.recognition, 'k--', 'Recognition req.');
xline(100*p.posErrorStd, 'k:', 'Design value');
xlabel('1\sigma position error (cm)'); ylabel('Worst run (%)'); grid on;
legend('Geometric coverage', 'Marker recognition', 'Location', 'southwest');
title('Sensitivity to position error');
subplot(1,3,2);
bar(categorical(compose('%.2f', accelValues)), accelSearchTime);
yline(req.searchTime, 'k--', 'Req. (TBD)');
xlabel('Commanded acceleration (m/s^2)'); ylabel('Search time (s)'); grid on;
title('Sensitivity to acceleration');
subplot(1,3,3);
plot(altitudes, markerPixelsVsAltitude, 'LineWidth', 1.5); hold on;
yline(p.minimumMarkerPixels, 'k--', 'Min pixels');
xline(p.altitude, 'k:', 'h_{max}');
xlabel('Altitude (m)'); ylabel('Pixels across marker'); grid on;
legend('720p', '1080p', '1440p', '4K', 'Location', 'northeast');
title(sprintf('Marker pixels vs altitude (%d^\\circ HFOV)', p.hFOV));
exportgraphics(gcf, fullfile(figDir, 'sensitivities.png'), 'Resolution', 150);

% Energy budget
newFigure([100 100 800 400]);
bar(categorical(budget.phases, budget.phases), budget.energy);
hold on; yline(req.energy, 'r--', sprintf('Allowed (%.0f%% reserve)', 100*e.reserveFraction));
ylabel('Energy (Wh)'); grid on;
title(sprintf('Mission energy budget: %.1f Wh of %.1f Wh allowed (inputs TBD)', budget.totalEnergy, req.energy));
exportgraphics(gcf, fullfile(figDir, 'energy_budget.png'), 'Resolution', 150);

%% Functions

function fp = cameraFootprint(altitude, p)
    fp.width = 2*altitude*tand(p.hFOV/2);    % cross-track
    fp.length = 2*altitude*tand(p.vFOV/2);   % along-track
    fp.gsd = fp.width/p.imageResolution(1);
end

% DR1 layout: tracks from edge to edge at footprint-width spacing
function waypoints = generateLawnmowerDR1(p, fp)
    trackSpacing = fp.width*(1 - p.overlap);
    xlines = 0:trackSpacing:p.areaWidth;
    if abs(xlines(end) - p.areaWidth) > 1e-10
        xlines = [xlines, p.areaWidth];
    end

    waypoints = [];
    for i = 1:numel(xlines)
        if mod(i,2) == 1
            ys = [0, p.areaLength];
        else
            ys = [p.areaLength, 0];
        end
        waypoints = [waypoints; xlines(i), ys(1); xlines(i), ys(2)];
    end
end

% DR2 layout: tracks inset so the footprint, shrunk by half a marker and the
% 3-sigma error margins, still reaches every possible marker centre
function waypoints = generateInsetLawnmower(p, fp, crossTrackMargin, alongTrackMargin)
    half = p.markerSize/2;
    halfWidth = fp.width/2 - half - crossTrackMargin;
    halfLength = fp.length/2 - half - alongTrackMargin;

    if halfWidth <= 0 || halfLength <= 0
        error("Error margins exceed the camera footprint")
    end

    span = p.areaWidth - 2*half;
    trackSpacing = 2*halfWidth*(1 - p.overlap);
    nTracks = max(1, ceil((span - 2*halfWidth)/trackSpacing) + 1);

    if nTracks == 1
        xlines = p.areaWidth/2;
    else
        xlines = linspace(half + halfWidth, p.areaWidth - half - halfWidth, nTracks);
    end

    yLow = min(half + halfLength, p.areaLength/2);
    yHigh = max(p.areaLength - half - halfLength, p.areaLength/2);

    waypoints = [];
    for i = 1:nTracks
        if mod(i,2) == 1
            ys = [yLow, yHigh];
        else
            ys = [yHigh, yLow];
        end
        waypoints = [waypoints; xlines(i), ys(1); xlines(i), ys(2)];
    end
end

% Commanded trajectory: each leg starts and ends at rest with a trapezoidal
% speed profile, followed by a settle hold at the waypoint
function traj = buildTrajectory(waypoints, speed, accel, cornerHold, dt)
    t = 0; pos = waypoints(1,:); vel = [0 0]; acc = [0 0];
    tNow = 0;

    for i = 1:size(waypoints,1)-1
        startPoint = waypoints(i,:);
        endPoint = waypoints(i+1,:);
        legLength = norm(endPoint - startPoint);
        if legLength < 1e-10, continue; end
        direction = (endPoint - startPoint)/legLength;

        [tau, s, sd, sdd] = trapezoidProfile(legLength, speed, accel, dt);
        tau = tau(2:end); s = s(2:end); sd = sd(2:end); sdd = sdd(2:end);

        t = [t; tNow + tau];
        pos = [pos; startPoint + s.*direction];
        vel = [vel; sd.*direction];
        acc = [acc; sdd.*direction];
        tNow = t(end);

        nHold = round(cornerHold/dt);
        if nHold > 0
            t = [t; tNow + (1:nHold)'*dt];
            pos = [pos; repmat(endPoint, nHold, 1)];
            vel = [vel; zeros(nHold, 2)];
            acc = [acc; zeros(nHold, 2)];
            tNow = t(end);
        end
    end

    traj.t = t; traj.pos = pos; traj.vel = vel; traj.acc = acc;
end

function [tau, s, sd, sdd] = trapezoidProfile(legLength, speed, accel, dt)
    tAccel = speed/accel;
    dAccel = 0.5*accel*tAccel^2;

    if 2*dAccel > legLength      % triangular: never reaches cruise speed
        tAccel = sqrt(legLength/accel);
        speed = accel*tAccel;
        tCruise = 0;
    else
        tCruise = (legLength - 2*dAccel)/speed;
    end

    T = 2*tAccel + tCruise;
    tau = (0:dt:T)';
    if T - tau(end) > 1e-9, tau = [tau; T]; end

    s = zeros(size(tau)); sd = s; sdd = s;
    phase1 = tau <= tAccel;
    phase2 = tau > tAccel & tau <= tAccel + tCruise;
    phase3 = tau > tAccel + tCruise;

    s(phase1) = 0.5*accel*tau(phase1).^2;
    sd(phase1) = accel*tau(phase1);
    sdd(phase1) = accel;

    s(phase2) = 0.5*accel*tAccel^2 + speed*(tau(phase2) - tAccel);
    sd(phase2) = speed;

    td = tau(phase3) - tAccel - tCruise;
    s(phase3) = 0.5*accel*tAccel^2 + speed*tCruise + speed*td - 0.5*accel*td.^2;
    sd(phase3) = speed - accel*td;
    sdd(phase3) = -accel;
end

% Leg time for a rest-to-rest trapezoidal profile (governing equation)
function T = legTime(legLength, speed, accel)
    if legLength >= speed^2/accel
        T = legLength/speed + speed/accel;
    else
        T = 2*sqrt(legLength/accel);
    end
end

function e = gaussMarkov(n, dim, sigma, tau, dt)
    phi = exp(-dt/tau);
    e = zeros(n, dim);
    e(1,:) = sigma*randn(1, dim);
    for k = 2:n
        e(k,:) = phi*e(k-1,:) + sigma*sqrt(1 - phi^2)*randn(1, dim);
    end
end

% One simulated search: sample frames, apply errors, score coverage and recognition
function res = simulateSearch(traj, p, X, Y, markerXY)
    dtFrame = 1/p.frameRate;
    tFrames = (0:dtFrame:traj.t(end))';
    nFrames = numel(tFrames);

    pos = interp1(traj.t, traj.pos, tFrames);
    vel = interp1(traj.t, traj.vel, tFrames);
    acc = interp1(traj.t, traj.acc, tFrames);

    posError = gaussMarkov(nFrames, 2, p.posErrorStd, p.posErrorTau, dtFrame);
    altError = gaussMarkov(nFrames, 1, p.altitudeErrorStd, p.posErrorTau, dtFrame);
    attJitter = gaussMarkov(nFrames, 2, deg2rad(p.attitudeJitterStd), p.attitudeTau, dtFrame);

    altitude = p.altitude + altError;

    % Accelerating tilts the body toward the acceleration, so a body-fixed
    % nadir camera looks opposite to it
    tiltLimit = deg2rad(p.tiltLimit);
    tilt = min(max(atan(acc/p.g), -tiltLimit), tiltLimit) + attJitter;
    centers = pos + posError - altitude.*tan(tilt);

    halfWidth = altitude*tand(p.hFOV/2);
    halfLength = altitude*tand(p.vFOV/2);
    gsd = 2*halfWidth/p.imageResolution(1);
    markerPixels = p.markerSize./gsd;
    blurPixels = vecnorm(vel, 2, 2)*p.exposure./gsd;
    goodFrame = markerPixels >= p.minimumMarkerPixels & blurPixels <= p.maxBlurPixels;

    countGrid = zeros(size(X));
    frameCount = zeros(size(markerXY,1), 1);
    recognitionTime = NaN(size(markerXY,1), 1);
    half = p.markerSize/2;

    for k = 1:nFrames
        inFootprint = abs(X - centers(k,1)) <= halfWidth(k) & abs(Y - centers(k,2)) <= halfLength(k);
        countGrid = countGrid + inFootprint;

        if goodFrame(k)
            inFrame = abs(markerXY(:,1) - centers(k,1)) <= halfWidth(k) - half & ...
                      abs(markerXY(:,2) - centers(k,2)) <= halfLength(k) - half;
            frameCount = frameCount + inFrame;
            newlyRecognized = isnan(recognitionTime) & frameCount >= p.minRecognitionFrames;
            recognitionTime(newlyRecognized) = tFrames(k);
        end
    end

    res.coverage = 100*mean(countGrid(:) > 0);
    res.recognition = 100*mean(frameCount >= p.minRecognitionFrames);
    res.worstRecognitionTime = max(recognitionTime);
    res.worstGSD = max(gsd);
    res.minPixels = min(markerPixels);
    res.maxBlur = max(blurPixels);
    res.countGrid = countGrid;
    res.centers = centers;
end

function mc = runMonteCarlo(traj, p, nRuns, X, Y, markerXY)
    mc.coverage = zeros(nRuns,1);
    mc.recognition = zeros(nRuns,1);
    mc.worstRecognitionTime = zeros(nRuns,1);
    mc.worstGSD = zeros(nRuns,1);
    mc.minPixels = zeros(nRuns,1);
    mc.maxBlur = zeros(nRuns,1);

    for r = 1:nRuns
        res = simulateSearch(traj, p, X, Y, markerXY);
        mc.coverage(r) = res.coverage;
        mc.recognition(r) = res.recognition;
        mc.worstRecognitionTime(r) = res.worstRecognitionTime;
        mc.worstGSD(r) = res.worstGSD;
        mc.minPixels(r) = res.minPixels;
        mc.maxBlur(r) = res.maxBlur;
        if r == 1
            mc.exampleCountGrid = res.countGrid;
            mc.exampleCenters = res.centers;
        end
    end
end

% Worst-case mission timeline and energy: full search before the marker is confirmed
function budget = missionBudget(traj, waypoints, p, m, e)
    transitTime = legTime(norm(waypoints(1,:) - m.launchPoint), p.vCruise, p.accel) + p.cornerHold;
    returnTime = legTime(norm(waypoints(end,:) - m.launchPoint), p.vCruise, p.accel) + p.cornerHold;
    landingTime = p.altitude/m.descentRate;

    budget.phases = {'Takeoff', 'Init hover', 'Transit', 'Area search', 'Marker hover', 'Comms hover', 'Return', 'Landing'};
    budget.time = [m.takeoffTime, m.initHover, transitTime, traj.t(end), m.markerHover, m.commsHover, returnTime, landingTime];

    % At <= 0.5 m/s, forward-flight power is approximately hover power
    climbPower = e.hoverPower + e.mass*p.g*(p.altitude/m.takeoffTime)/e.climbEfficiency;
    power = e.hoverPower*ones(size(budget.time));
    power(1) = climbPower;

    budget.energy = power.*budget.time/3600; % [Wh]
    budget.totalTime = sum(budget.time);
    budget.totalEnergy = sum(budget.energy);
end

function printRow(name, requirement, value, pass)
    if pass, status = 'PASS'; else, status = 'FAIL'; end
    if strcmp(requirement, '-'), status = '-'; end
    fprintf('%-34s %-14s %-14.2f %-8s\n', name, requirement, value, status);
end

% Light-theme figure so plots read correctly on slides
function fig = newFigure(position)
    fig = figure('Position', position, 'Color', 'w');
    try
        theme(fig, 'light');
    catch
    end
end

function plotPath(waypoints, p, titleText)
    rectangle('Position', [0 0 p.areaWidth p.areaLength], 'FaceColor', [0.96 0.96 0.96], 'EdgeColor', 'k'); hold on;
    plot(waypoints(:,1), waypoints(:,2), '-o', 'LineWidth', 1.5, 'MarkerSize', 4);
    plot(waypoints(1,1), waypoints(1,2), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 8);
    plot(waypoints(end,1), waypoints(end,2), 'rs', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
    axis equal; xlim([-0.3 p.areaWidth+0.3]); ylim([-0.3 p.areaLength+0.3]); grid on;
    xlabel('X Position (m)'); ylabel('Y Position (m)');
    title(titleText);
end
