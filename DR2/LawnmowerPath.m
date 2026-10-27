%% DR2 inset lawnmower layout (single plot)
clc, clear, close all;

scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir), scriptDir = pwd; end
figDir = fullfile(scriptDir, 'figures');
if ~exist(figDir, 'dir'), mkdir(figDir); end

%% Search area and flight envelope
p.areaWidth = 5;              % [m]
p.areaLength = 5;             % [m]
p.altitude = 2.0;             % [m]
p.vCruise = 0.40;             % [m/s]
p.g = 9.81;                   % [m/s^2]

%% Path-following assumptions
p.accel = 0.50;               % [m/s^2]
p.cornerHold = 1.0;           % [s]
p.overlap = 0.20;             % [-]
p.posErrorStd = 0.05;         % [m]
p.attitudeJitterStd = 1.0;    % [deg]
p.dt = 0.01;                  % [s]

%% Camera
p.imageResolution = [3264, 2448]; % [px]
p.hFOV = 60;                  % [deg]
p.vFOV = 2*atand((p.imageResolution(2)/p.imageResolution(1))*tand(p.hFOV/2));

%% Marker
p.markerSize = 0.08;          % [m]

launchPoint = [-2, 0];        % [m]

%% Footprint and error margins
fp = cameraFootprint(p.altitude, p);
crossTrackMargin = 3*p.posErrorStd + p.altitude*tand(3*p.attitudeJitterStd);
alongTrackMargin = crossTrackMargin + p.altitude*(p.accel/p.g);

%% DR2 path
wpDR2 = generateInsetLawnmower(p, fp, crossTrackMargin, alongTrackMargin);

entryPoint = min(max(launchPoint, [0 0]), [p.areaWidth p.areaLength]);
wpDR2 = [entryPoint; wpDR2];

lastPoint = wpDR2(end,:);
edges = [0, p.areaWidth, 0, p.areaLength];
[~, k] = min(abs(lastPoint([1 1 2 2]) - edges));
exitPoint = lastPoint;
exitPoint(ceil(k/2)) = edges(k);
wpDR2 = [wpDR2; exitPoint];

trajDR2 = buildTrajectory(wpDR2, p.vCruise, p.accel, p.cornerHold, p.dt);
pathLengthDR2 = sum(vecnorm(diff(wpDR2), 2, 2));

%% Plot
newFigure([100 100 520 480]);
plotPath(wpDR2, p, sprintf('DR2 inset layout: %.1f m, %.1f s', pathLengthDR2, trajDR2.t(end)));
exportgraphics(gcf, fullfile(figDir, 'dr2_layout.png'), 'Resolution', 150);

%% Functions

function fp = cameraFootprint(altitude, p)
    fp.width = 2*altitude*tand(p.hFOV/2);
    fp.length = 2*altitude*tand(p.vFOV/2);
end

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

function traj = buildTrajectory(waypoints, speed, accel, cornerHold, dt)
    t = 0; tNow = 0;

    for i = 1:size(waypoints,1)-1
        legLength = norm(waypoints(i+1,:) - waypoints(i,:));
        if legLength < 1e-10, continue; end

        [tau, ~, ~, ~] = trapezoidProfile(legLength, speed, accel, dt);
        tau = tau(2:end);
        t = [t; tNow + tau];
        tNow = t(end);

        nHold = round(cornerHold/dt);
        if nHold > 0
            t = [t; tNow + (1:nHold)'*dt];
            tNow = t(end);
        end
    end

    traj.t = t;
end

function [tau, s, sd, sdd] = trapezoidProfile(legLength, speed, accel, dt)
    tAccel = speed/accel;
    dAccel = 0.5*accel*tAccel^2;

    if 2*dAccel > legLength
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
end

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