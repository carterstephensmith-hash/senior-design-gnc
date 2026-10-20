clear; clc; close all;

%% Simulation timing
imuRate = 200;       % IMU measurements used per second [Hz]
rtkRate = 20;        % RTK measurements used per second [Hz]

dt = 1/imuRate;      % Time between simulation steps [s]
rtkPeriod = imuRate/rtkRate;

%% Lawnmower path settings
areaWidth = 5;       % Search-area width [m]
areaLength = 5;      % Search-area length [m]
altitude = 2;        % Nominal search altitude [m]

hFOV = 60;          % Assumed horizontal camera field of view [deg]
overlap = 0.20;     % Desired overlap between neighboring passes

maxSpeed = 0.5;     % Maximum aircraft speed [m/s]
maxAcceleration = 0.25; % Assumed acceleration/deceleration limit [m/s^2]

%% Monte Carlo settings

nTrials = 100;

% Column order for the stored results
caseNames = {'Direct, normal', 'Kalman, normal',  'Direct, RTK outage', 'Kalman, RTK outage'};

nCases = numel(caseNames);

% Each row represents one simulated flight
positionRMSE_mc = zeros(nTrials, nCases);
maxPositionError_mc = zeros(nTrials, nCases);

%% Calculate spacing between search passes
footprintWidth = 2*altitude*tand(hFOV/2);
trackSpacing = footprintWidth*(1-overlap);

%% Create lawnmower waypoints
waypoints = generateLawnmowerPath(areaWidth, areaLength, trackSpacing, altitude);

%% Generate the true aircraft motion
[t, p_true, v_true, a_true] = sampleLawnmowerPath(waypoints, maxSpeed, maxAcceleration, dt);

N = length(t);

%% Inspect the true aircraft motion
plotTrueMotion(t, p_true, v_true, a_true, maxSpeed, maxAcceleration);

%% Sample the true motion at RTK measurement times
[t_rtk, p_rtk_ideal, v_rtk_ideal, rtkIndices] = sampleRtkTruth(t, p_true, v_true, rtkPeriod);

%% RTK outage scenario

outageStart = 15; % seconds
outageEnd = 20;   % seconds

% Normal-operation case: every RTK measurement is available
rtkAvailable_normal = true(size(t_rtk));

% Outage case: measurements are unavailable from 15 up to 20 seconds
rtkAvailable_dropout = ~(t_rtk >= outageStart & t_rtk < outageEnd);

%% RTK measurement uncertainty
% Assumptions: zero-mean Gaussian errors, independent between axes.

% Horizontal accuracy: nominal 50% probability radius [m]
rtkHorizontalCEP = 0.01;

% Standard deviation of each horizontal component [m]
% Derive this equation by integrating the two dimensional Gaussian distribution over a circle of radius C
rtkHorizontalStd = rtkHorizontalCEP/sqrt(2*log(2));

% Vertical accuracy: nominal median absolute error [m]
rtkVerticalMedian = 0.01;

% Standard-normal value containing 75% of probability to its left
normalQuantile75 = sqrt(2)*erfinv(0.5);

% Vertical standard deviation [m]
rtkVerticalStd = rtkVerticalMedian/normalQuantile75;

% Position standard deviations ordered [East; North; Up]
rtkPositionStd = [rtkHorizontalStd;
                  rtkHorizontalStd;
                  rtkVerticalStd];

% Provisional velocity standard deviation for each axis [m/s]
% Retained from the original script; statistical meaning needs justification.
rtkVelocityStd = 0.05;

%% Fixed IMU accelerometer-noise settings

% Noise densities ordered [sensor X; sensor Y; sensor Z]
accelNoiseDensity_ug = [65; 65; 70]; % micro-g/sqrt(Hz)

% Convert to m/s^2/sqrt(Hz)
accelNoiseDensity = accelNoiseDensity_ug * 1e-6 * 9.81;

% Assumed equivalent noise bandwidth
accelNoiseBandwidth = 100; % Hz

accelNoiseStd = accelNoiseDensity * sqrt(accelNoiseBandwidth);

%% Fixed accelerometer-bias settings

% Research-based assumptions; validate against hardware data.
accelBiasStd = 0.002; % m/s^2
biasTau = 40;        % seconds

%% Ideal specific force and true orientation
% These depend on the prescribed trajectory, not the random trial.

gravityMission = [0; 0; -9.81]; % m/s^2, ENU frame

f_true_mission = a_true - gravityMission;

specificForceMagnitude = sqrt(sum(f_true_mission.^2, 1));
thrustDirectionMission = f_true_mission ./ specificForceMagnitude;
bodyZMission = -thrustDirectionMission;

R_bodyToMission = generateTrueAttitude(bodyZMission);

f_true_body = zeros(3, N);

for k = 1:N
    f_true_body(:, k) = R_bodyToMission(:, :, k)' * f_true_mission(:, k);
end

%% Fixed attitude-uncertainty settings

% Exploratory assumptions, not verified Pixhawk specifications.
attitudeErrorStdDeg = [1; 1; 3]; % degrees
attitudeErrorStd = deg2rad(attitudeErrorStdDeg);

%% Requirement and filter uncertainty settings

positionRequirement = 0.10; % 3D position-error limit, m

initialStateStd = [rtkPositionStd;
                   rtkVelocityStd * ones(3, 1);
                   accelBiasStd * ones(3, 1);
                   attitudeErrorStd];

P0_kf = diag(initialStateStd.^2);

[H, R_rtk] = buildRtkMeasurementModel(rtkPositionStd, rtkVelocityStd);

%% Monte Carlo trials

% One row per case, one column per simulation time
sumSquaredPositionError_mc = zeros(nCases, N);

% Store East and North coordinates for every trial
paths_direct_mc = zeros(2, N, nTrials);
paths_kf_mc = zeros(2, N, nTrials);

paths_direct_rtkOut_mc = zeros(2, N, nTrials);
paths_kf_rtkOut_mc = zeros(2, N, nTrials);

for trial = 1:nTrials

    % Reproducible, different random errors for each trial
    rng(trial);

    %% Generate one shared set of sensor errors for this trial

    [p_rtk, v_rtk] = generateRtkMeasurements(p_rtk_ideal, v_rtk_ideal, rtkPositionStd, rtkVelocityStd);

    b_accel = generateAccelBias(N, dt, accelBiasStd, biasTau);

    f_imu_body = generateAccelerometerMeasurements(f_true_body, b_accel, accelNoiseStd);

    % One persistent attitude offset per trial
    attitudeErrorAngles = attitudeErrorStd .* randn(3, 1);

    R_est_bodyToMission = generateEstimatedAttitude(R_bodyToMission, attitudeErrorAngles);

    a_imu_mission = convertAccelToMission(f_imu_body, R_est_bodyToMission, gravityMission);

    %% Initialize the filter using this trial's first RTK measurement

    x0_kf = [p_rtk(:, 1);
             v_rtk(:, 1);
             zeros(3, 1);
             zeros(3, 1)];

    %% Case 1: Direct RTK + IMU, normal operation

    [p_direct, v_direct] = runDirectRtkImu(a_imu_mission, p_rtk, v_rtk, rtkIndices, dt, rtkAvailable_normal);

    [positionError_direct, positionRMSE_direct, maxPositionError_direct] = evaluatePositionError(p_direct, p_true);

    %% Case 2: Kalman filter, normal operation

    [p_kf, v_kf, b_kf] = runKalmanFilter( a_imu_mission, f_imu_body, R_est_bodyToMission, p_rtk, v_rtk, rtkIndices, dt, ...
        x0_kf, P0_kf, H, R_rtk, accelNoiseStd, accelBiasStd, biasTau, rtkAvailable_normal);

    [positionError_kf, positionRMSE_kf, maxPositionError_kf] = evaluatePositionError(p_kf, p_true);

    %% Case 3: Direct RTK + IMU, RTK outage

    [p_direct_rtkOut, v_direct_rtkOut] = runDirectRtkImu(a_imu_mission, p_rtk, v_rtk, rtkIndices, dt, rtkAvailable_dropout);

    [positionError_direct_rtkOut, rmseDirectAll, maxDirectAll] = evaluatePositionError(p_direct_rtkOut, p_true);

    %% Case 4: Kalman filter, RTK outage

    [p_kf_rtkOut, v_kf_rtkOut, b_kf_rtkOut] = runKalmanFilter(a_imu_mission, f_imu_body, R_est_bodyToMission, ...
        p_rtk, v_rtk, rtkIndices, dt, x0_kf, P0_kf, H, R_rtk, accelNoiseStd, accelBiasStd, biasTau, rtkAvailable_dropout);

    [positionError_kf_rtkOut, rmseKfAll, maxKfAll] = evaluatePositionError(p_kf_rtkOut, p_true);

    %% Store whole-flight results for this trial
    % Columns:
    % 1 = Direct, normal
    % 2 = Kalman, normal
    % 3 = Direct, RTK outage
    % 4 = Kalman, RTK outage

    positionRMSE_mc(trial, :) = [positionRMSE_direct, positionRMSE_kf, rmseDirectAll, rmseKfAll];

    maxPositionError_mc(trial, :) = [maxPositionError_direct, maxPositionError_kf, maxDirectAll, maxKfAll];

    trialPositionErrors = [positionError_direct;
                       positionError_kf;
                       positionError_direct_rtkOut;
                       positionError_kf_rtkOut];

    sumSquaredPositionError_mc = sumSquaredPositionError_mc + trialPositionErrors.^2;

    paths_direct_mc(:, :, trial) = p_direct(1:2, :);
    paths_kf_mc(:, :, trial) = p_kf(1:2, :);

    paths_direct_rtkOut_mc(:, :, trial) = p_direct_rtkOut(1:2, :);
    paths_kf_rtkOut_mc(:, :, trial) = p_kf_rtkOut(1:2, :);

    %% Display detailed results for trial 1 only
    % setting trial == 0 turns this if statement off for now
    if trial == 0

        plotRtkErrors( t_rtk, p_rtk, p_rtk_ideal, v_rtk, v_rtk_ideal);

        % Direct-estimator position error
        figure;
        plot(t, positionError_direct, 'LineWidth', 1.2);
        hold on;
        yline(positionRequirement, '--r', '10 cm requirement');

        xlabel('Time (s)');
        ylabel('3D position error (m)');
        title('Direct RTK + IMU Position Error — Trial 1');
        grid on;

        % Normal-operation comparison
        figure;
        plot(t, positionError_direct, 'LineWidth', 1.2);
        hold on;
        plot(t, positionError_kf, 'LineWidth', 1.2);
        yline(positionRequirement, '--r');

        xlabel('Time (s)');
        ylabel('3D position error (m)');
        title('Localization Position Error — Trial 1');
        legend('Direct RTK + IMU', 'Kalman filter', '10 cm requirement', 'Location', 'best');
        grid on;

        % RTK-outage comparison
        figure;
        plot(t, positionError_direct_rtkOut, 'LineWidth', 1.2);
        hold on;
        plot(t, positionError_kf_rtkOut, 'LineWidth', 1.2);
        yline(positionRequirement, '--r');

        xline(outageStart, ':k', 'RTK lost', 'HandleVisibility', 'off');
        xline(outageEnd, ':k', 'RTK restored', 'HandleVisibility', 'off');

        xlabel('Time (s)');
        ylabel('3D position error (m)');
        title('Position Error During an RTK Outage — Trial 1');
        legend('Direct RTK + IMU', 'Kalman filter', '10 cm requirement', 'Location', 'best');
        grid on;

        % Optional zoom; larger errors would be outside this view:
        % ylim([0 0.5]);

        % Whole-flight metrics for trial 1
        fprintf('\nTrial 1: normal operation\n');
        fprintf('Direct: RMSE = %.4f m, maximum = %.4f m\n', positionRMSE_direct, maxPositionError_direct);
        fprintf('Kalman: RMSE = %.4f m, maximum = %.4f m\n', positionRMSE_kf, maxPositionError_kf);

        fprintf('\nTrial 1: RTK outage from %.1f to %.1f s\n', outageStart, outageEnd);
        fprintf('Direct: RMSE = %.4f m, maximum = %.4f m\n', rmseDirectAll, maxDirectAll);
        fprintf('Kalman: RMSE = %.4f m, maximum = %.4f m\n', rmseKfAll, maxKfAll);

    end

    % Progress update
    if mod(trial, 10) == 0 || trial == nTrials
        fprintf('Completed trial %d of %d.\n', trial, nTrials);
    end

end

%% Monte Carlo results

% RMS position error across trials at each time
rmsPositionErrorTime_mc = sqrt(sumSquaredPositionError_mc / nTrials);

% Whole-flight RMSE pooled across all trials and time samples
pooledPositionRMSE_mc = sqrt(mean(positionRMSE_mc.^2, 1));

% Maximum-error summaries
meanTrialMaximum_mc = mean(maxPositionError_mc, 1);
worstObservedError_mc = max(maxPositionError_mc, [], 1);

% Percentage of trials remaining within the requirement for the entire flight
percentTrialsWithinRequirement_mc = 100 * mean(maxPositionError_mc <= positionRequirement, 1);

% Monte Carlo normal-operation plot

figure;
plot(t, rmsPositionErrorTime_mc(1, :), 'LineWidth', 1.5);
hold on;
plot(t, rmsPositionErrorTime_mc(2, :), 'LineWidth', 1.5);
yline(positionRequirement, '--r');

xlabel('Time (s)');
ylabel('RMS 3D position error across trials (m)');
title(sprintf('Normal Operation — %d Monte Carlo Trials', nTrials));
legend('Direct RTK + IMU', 'Kalman filter', '10 cm reference', 'Location', 'best');
grid on;

% Monte Carlo RTK-outage plot

figure;
plot(t, rmsPositionErrorTime_mc(3, :), 'LineWidth', 1.5);
hold on;
plot(t, rmsPositionErrorTime_mc(4, :), 'LineWidth', 1.5);
yline(positionRequirement, '--r');

xline(outageStart, ':k', 'RTK lost', 'HandleVisibility', 'off');
xline(outageEnd, ':k', 'RTK restored', 'HandleVisibility', 'off');

xlabel('Time (s)');
ylabel('RMS 3D position error across trials (m)');
title(sprintf('RTK Outage — %d Monte Carlo Trials', nTrials));
legend('Direct RTK + IMU', 'Kalman filter', '10 cm reference', 'Location', 'best');
grid on;

% Print aggregate metrics

fprintf('\nMonte Carlo results: %d trials\n', nTrials);

for caseIndex = 1:nCases

    fprintf('\n%s\n', caseNames{caseIndex});
    fprintf('  Pooled position RMSE: %.4f m\n', pooledPositionRMSE_mc(caseIndex));
    fprintf('  Mean trial maximum error: %.4f m\n', meanTrialMaximum_mc(caseIndex));
    fprintf('  Worst observed error: %.4f m\n', worstObservedError_mc(caseIndex));
    fprintf('  Trials entirely within 10 cm: %.1f%%\n', percentTrialsWithinRequirement_mc(caseIndex));

end

plotMonteCarloPaths(waypoints, p_true, paths_direct_mc, paths_kf_mc, positionRequirement, sprintf('Normal Operation — %d Monte Carlo Trials', nTrials));

plotMonteCarloPaths(waypoints, p_true, paths_direct_rtkOut_mc, paths_kf_rtkOut_mc, positionRequirement, ...
    sprintf('RTK Outage %.0f–%.0f s — %d Monte Carlo Trials', outageStart, outageEnd, nTrials));

function plotTrueMotion(t, p_true, v_true, a_true, maxSpeed, maxAcceleration)

    % Magnitudes of the 3D velocity and acceleration vectors
    trueSpeed = sqrt(sum(v_true.^2,1));
    trueAccelMagnitude = sqrt(sum(a_true.^2,1));

    % Plot the flight path
    figure;
    plot3(p_true(1,:), p_true(2,:), p_true(3,:), 'b-');
    grid on;
    axis equal;
    xlabel('East [m]');
    ylabel('North [m]');
    zlabel('Altitude [m]');
    title('True aircraft trajectory');

    % Plot speed and acceleration
    figure;

    subplot(2,1,1);
    plot(t, trueSpeed, 'b-');
    yline(maxSpeed, 'r--', 'Speed limit');
    grid on;
    xlabel('Time [s]');
    ylabel('Speed [m/s]');
    title('True aircraft speed');

    subplot(2,1,2);
    plot(t, trueAccelMagnitude, 'b-');
    yline(maxAcceleration, 'r--', 'Acceleration limit');
    grid on;
    xlabel('Time [s]');
    ylabel('Acceleration magnitude [m/s^2]');
    title('True aircraft acceleration magnitude');

    % Print useful checks
    fprintf('Simulation duration: %.3f s\n', t(end));
    fprintf('Maximum speed: %.3f m/s\n', max(trueSpeed));
    fprintf('Maximum acceleration: %.3f m/s^2\n', ...
        max(trueAccelMagnitude));
end

function [t_rtk, p_rtk_ideal, v_rtk_ideal, rtkIndices] = sampleRtkTruth(t, p_true, v_true, rtkPeriod)

    % Indices corresponding to RTK measurement times
    rtkIndices = 1:rtkPeriod:length(t);

    % Select the time, position, and velocity at those indices
    t_rtk = t(rtkIndices);
    p_rtk_ideal = p_true(:,rtkIndices);
    v_rtk_ideal = v_true(:,rtkIndices);
end

function [p_rtk, v_rtk] = generateRtkMeasurements(p_rtk_ideal, v_rtk_ideal, rtkPositionStd, rtkVelocityStd)

    % Generate independent Gaussian measurement errors
    positionNoise = rtkPositionStd .* randn(size(p_rtk_ideal));
    velocityNoise = rtkVelocityStd .* randn(size(v_rtk_ideal));

    % Add errors to the ideal measurements
    p_rtk = p_rtk_ideal + positionNoise;
    v_rtk = v_rtk_ideal + velocityNoise;
end

function plotRtkErrors(t_rtk, p_rtk, p_rtk_ideal, v_rtk, v_rtk_ideal)

    % Measurement error = reported value minus true value
    positionError = p_rtk - p_rtk_ideal;
    velocityError = v_rtk - v_rtk_ideal;

    figure;

    subplot(2,1,1);
    plot(t_rtk, positionError');
    grid on;
    xlabel('Time [s]');
    ylabel('Position error [m]');
    title('Simulated RTK position errors');
    legend('East', 'North', 'Up');

    subplot(2,1,2);
    plot(t_rtk, velocityError');
    grid on;
    xlabel('Time [s]');
    ylabel('Velocity error [m/s]');
    title('Simulated RTK velocity errors');
    legend('East', 'North', 'Up');

    % Calculate the sample standard deviation along each axis
    fprintf('Position error standard deviations [East; North; Up]:\n');
    disp(std(positionError, 0, 2));

    fprintf('Velocity error standard deviations [East; North; Up]:\n');
    disp(std(velocityError, 0, 2));
end

% First order Guass Markov bias Model with an assumed magnitude and correlation time
function b_accel = generateAccelBias(N, dt, accelBiasStd, biasTau)

    % Fraction of the previous bias retained at each time step
    alpha = exp(-dt / biasTau);

    % Standard deviation of the new random contribution
    biasStepStd = accelBiasStd * sqrt(1 - alpha^2);

    % One row per accelerometer axis, one column per time step
    b_accel = zeros(3, N);

    % Random initial bias
    b_accel(:, 1) = accelBiasStd * randn(3, 1);

    % Generate the remaining bias history
    for k = 2:N
        b_accel(:, k) = alpha * b_accel(:, k-1) + biasStepStd * randn(3, 1);
    end
end

function R_bodyToMission = generateTrueAttitude(bodyZMission)

    N = size(bodyZMission, 2);
    R_bodyToMission = zeros(3, 3, N);

    % With the nose pointing east, horizontal right is south
    referenceRightMission = [0; -1; 0];

    for k = 1:N

        zBody = bodyZMission(:, k);

        % Construct the forward axis
        xBody = cross(referenceRightMission, zBody);
        xBody = xBody / norm(xBody);

        % Construct the actual right axis
        yBody = cross(zBody, xBody);

        % Columns are the body axes expressed in mission coordinates
        R_bodyToMission(:, :, k) = [xBody, yBody, zBody];

    end

end

function f_imu_body = generateAccelerometerMeasurements(f_true_body, b_accel, accelNoiseStd)

    % Generate independent random measurement noise
    accelNoise = accelNoiseStd .* randn(size(f_true_body));

    % Combine ideal measurement, slowly varying bias, and noise
    f_imu_body = f_true_body + b_accel + accelNoise;

end

function R_est_bodyToMission = generateEstimatedAttitude(R_bodyToMission, attitudeErrorAngles)

    % Error angles are already in radians
    ex = attitudeErrorAngles(1);
    ey = attitudeErrorAngles(2);
    ez = attitudeErrorAngles(3);

    % Rotation about the x axis
    Rx = [1,       0,        0;
          0, cos(ex), -sin(ex);
          0, sin(ex),  cos(ex)];

    % Rotation about the y axis
    Ry = [ cos(ey), 0, sin(ey);
                 0, 1,       0;
          -sin(ey), 0, cos(ey)];

    % Rotation about the z axis
    Rz = [cos(ez), -sin(ez), 0;
          sin(ez),  cos(ez), 0;
                0,        0, 1];

    % Combine the three error rotations
    R_error = Rz * Ry * Rx;

    N = size(R_bodyToMission, 3);
    R_est_bodyToMission = zeros(3, 3, N);

    for k = 1:N
        R_est_bodyToMission(:, :, k) = R_bodyToMission(:, :, k) * R_error;
    end

end

function a_imu_mission = convertAccelToMission(f_imu_body, R_est_bodyToMission, gravityMission)

    N = size(f_imu_body, 2);
    a_imu_mission = zeros(3, N);

    for k = 1:N

        % Rotate measured specific force into the mission frame
        f_mission_est = R_est_bodyToMission(:, :, k) * f_imu_body(:, k);

        % Recover motion acceleration by accounting for gravity
        a_imu_mission(:, k) = f_mission_est + gravityMission;

    end

end

function [H, R_rtk] = buildRtkMeasurementModel(rtkPositionStd, rtkVelocityStd)

    % RTK observes position and velocity, not bias directly
    H = [eye(6), zeros(6, 6)];

    measurementStd = [rtkPositionStd;
                      rtkVelocityStd * ones(3, 1)];

    R_rtk = diag(measurementStd.^2);

end

