clear; clc; close all;

%% Marker geometry
markerSideLength = 0.203; % Outer black-square side length [m], provisional

%% Nominal camera settings: original Arducam B029201
imageWidth = 3264;        % Initial modeled image width [pixels]
imageHeight = 2448;       % Initial modeled image height [pixels]
horizontalFovDeg = 60;    % Manufacturer's nominal horizontal FOV [deg]

%% Ideal-camera assumptions for the geometry check
fx = (imageWidth/2) / tand(horizontalFovDeg/2); % [pixels]
fy = fx;                                     % Equal image scaling

% MATLAB image convention: first pixel center is at coordinate 1
cx = (imageWidth + 1)/2;   % Principal-point horizontal coordinate [pixels]
cy = (imageHeight + 1)/2;  % Principal-point vertical coordinate [pixels]

%% True marker geometry in the mission frame
markerCenterMission = [2.5; 2.5; 0]; % [East; North; Up], m

halfSide = markerSideLength/2;

% Each column describes one corner relative to the marker center.
% Components are along mission East, North, and Up.
markerCornerOffsetsMission = halfSide * [-1,  1,  1, -1;
                                         -1, -1,  1,  1;
                                          0,  0,  0,  0];

markerCornersMission = markerCenterMission + markerCornerOffsetsMission;

%% Assumed nominal camera mounting orientation
% Columns are the camera axes expressed in body coordinates.
% Xc image right -> Yb aircraft right
% Yc image down -> -Xb aircraft backward
% Zc outward through lens -> Zb aircraft down
R_cameraToBody = [0, -1, 0;
                 1,  0, 0;
                 0,  0, 1];

%% Complete ideal-geometry verification

% Known physical marker corners, relative to the marker center [m].
% Row order must match the projected corner order.
markerCornersLocal = halfSide * [-1, -1;
                                  1, -1;
                                  1,  1;
                                 -1,  1];

intrinsics = cameraIntrinsics([fx, fy], [cx, cy], [imageHeight, imageWidth]);

bodyPositionMission = [2; 2; 2]; % [East; North; Up], m

% Columns: roll [deg], pitch [deg], yaw from east [deg],
%          camera offset along body-down [m].
testCases = [0,  0,  0, 0.00;
             0,  0,  0, 0.10;
             8, -5, 20, 0.10];

caseNames = {'Level, zero offset';
             'Level, camera offset';
             'Tilted, camera offset'};

% Columns: camera-frame, body-frame, and mission-frame 3D errors [m].
geometryErrors = zeros(size(testCases, 1), 3);

figure;
hold on;

for caseIndex = 1:size(testCases, 1)

    %% 1. Define the true aircraft pose and camera mounting
    R_bodyToMission = bodyToMissionFromAngles(testCases(caseIndex, 1), testCases(caseIndex, 2), testCases(caseIndex, 3));

    cameraPositionBody = [0; 0; testCases(caseIndex, 4)];

    %% 2. Generate the true camera-frame corner positions
    markerCornersBody = R_bodyToMission' * (markerCornersMission - bodyPositionMission);

    markerCornersCamera = R_cameraToBody' * (markerCornersBody - cameraPositionBody);

    Xc = markerCornersCamera(1,:);
    Yc = markerCornersCamera(2,:);
    Zc = markerCornersCamera(3,:);

    assert(all(Zc > 0), 'A marker corner is behind the camera in case %d.', caseIndex);

    %% 3. Generate perfect image measurements
    u = fx * (Xc ./ Zc) + cx;
    v = fy * (Yc ./ Zc) + cy;

    insideImage = (u >= 0.5 & u <= imageWidth + 0.5) & (v >= 0.5 & v <= imageHeight + 0.5);

    assert(all(insideImage), 'A marker corner is outside the image in case %d.', caseIndex);

    markerPixels = [u; v];

    %% 4. Recover position using image measurements and known marker size
    markerToCamera = estimateExtrinsics(markerPixels', markerCornersLocal, intrinsics);

    markerCenterCamera_est = markerToCamera.Translation';

    %% 5. Convert the recovered position into body and mission coordinates
    markerCenterBody_est = cameraPositionBody + R_cameraToBody * markerCenterCamera_est;

    markerCenterMission_est = bodyPositionMission + R_bodyToMission * markerCenterBody_est;

    %% 6. Evaluate against truth, separately from the estimator
    markerCenterBody_true = R_bodyToMission' * (markerCenterMission - bodyPositionMission);

    markerCenterCamera_true = R_cameraToBody' * (markerCenterBody_true - cameraPositionBody);

    geometryErrors(caseIndex,:) = [norm(markerCenterCamera_est - markerCenterCamera_true), norm(markerCenterBody_est - markerCenterBody_true), ...
        norm(markerCenterMission_est - markerCenterMission)];

    %% Plot this case's ideal image measurements
    closedOrder = [1:4, 1];

    plot(u(closedOrder), v(closedOrder), '-o', 'LineWidth', 1.5);
end

axis equal;
xlim([0.5, imageWidth + 0.5]);
ylim([0.5, imageHeight + 0.5]);
set(gca, 'YDir', 'reverse');
grid on;
xlabel('u [pixels]');
ylabel('v [pixels]');
title('Ideal marker projections for three verification cases');
legend(caseNames, 'Location', 'best');

geometryResults = table(caseNames, geometryErrors(:,1), geometryErrors(:,2), geometryErrors(:,3), 'VariableNames', {'Case', 'CameraError_m', ...
                     'BodyError_m', 'MissionError_m'});

disp(geometryResults);

%% Generate operating-condition geometries

% Initial study points, not verified flight-envelope limits.
studyAltitudes_m = [1, 1.5, 2];

% Columns: roll, pitch, yaw [deg].
studyAttitudes_deg = [0, 0, 0;
                     5, 0, 0;
                     0, 5, 0];

studyAttitudeNames = {'Level', '5 deg roll', '5 deg pitch'};

% Marker-center locations: center, near-center, intermediate, and outer.
studyImageFractions = [ ...
    0.5, 0.5;
    0.4, 0.4;
    0.6, 0.4;
    0.6, 0.6;
    0.4, 0.6;
    0.3, 0.3;
    0.7, 0.3;
    0.7, 0.7;
    0.3, 0.7;
    0.2, 0.2;
    0.8, 0.2;
    0.8, 0.8;
    0.2, 0.8];

studyLocationNames = { ...
    'Center', ...
    'Near upper left', ...
    'Near upper right', ...
    'Near lower right', ...
    'Near lower left', ...
    'Middle upper left', ...
    'Middle upper right', ...
    'Middle lower right', ...
    'Middle lower left', ...
    'Outer upper left', ...
    'Outer upper right', ...
    'Outer lower right', ...
    'Outer lower left'};

% Reuse the mounting offset from the existing camera-offset test case.
% This remains a provisional installation geometry.
studyCameraPositionBody = [0; 0; testCases(2,4)];

nAltitudes = numel(studyAltitudes_m);
nAttitudes = size(studyAttitudes_deg, 1);
nLocations = size(studyImageFractions, 1);

operatingCases = struct([]);
studyCaseIndex = 0;

figure;
tiledlayout(nAttitudes, nAltitudes, ...
    'TileSpacing', 'compact', 'Padding', 'compact');

for attitudeIndex = 1:nAttitudes

    angles = studyAttitudes_deg(attitudeIndex,:);

    studyR = bodyToMissionFromAngles(angles(1), angles(2), angles(3));

    % Camera mounting offset expressed in mission coordinates [m].
    cameraOffsetMission = studyR * studyCameraPositionBody;

    for altitudeIndex = 1:nAltitudes

        altitude_m = studyAltitudes_m(altitudeIndex);

        ax = nexttile;
        hold(ax, 'on');

        for locationIndex = 1:nLocations

            studyCaseIndex = studyCaseIndex + 1;

            % Desired pixel coordinates of the marker center.
            targetPixel = [0.5 + imageWidth * studyImageFractions(locationIndex,1);
                           0.5 + imageHeight * studyImageFractions(locationIndex,2)];

            % Camera viewing ray corresponding to that pixel.
            % Its components are dimensionless; it is not a unit vector.
            rayCamera = [(targetPixel(1) - cx)/fx;
                         (targetPixel(2) - cy)/fy;
                         1];

            rayMission = studyR * R_cameraToBody * rayCamera;

            % Select the aircraft position that puts the marker on this ray
            % while maintaining the specified body-origin altitude.
            rayScale = (markerCenterMission(3) - altitude_m - cameraOffsetMission(3)) / rayMission(3);

            assert(rayScale > 0, 'The selected geometry points away from the marker.');

            studyBodyPosition = markerCenterMission - cameraOffsetMission - rayScale * rayMission;

            % Forward-project the true marker corners.
            cornersBody = studyR' * (markerCornersMission - studyBodyPosition);

            cornersCamera = R_cameraToBody' * (cornersBody - studyCameraPositionBody);

            studyPixels = [fx * cornersCamera(1,:) ./ cornersCamera(3,:) + cx;
                           fy * cornersCamera(2,:) ./ cornersCamera(3,:) + cy];

            insideImage = studyPixels(1,:) >= 0.5 & studyPixels(1,:) <= imageWidth + 0.5 & studyPixels(2,:) >= 0.5 & studyPixels(2,:) <= imageHeight + 0.5;

            % Store the condition for the later Monte Carlo study.
            operatingCases(studyCaseIndex).name = sprintf('h=%.1f m, %s, %s', altitude_m, studyAttitudeNames{attitudeIndex}, studyLocationNames{locationIndex});

            operatingCases(studyCaseIndex).altitude_m = altitude_m;
            operatingCases(studyCaseIndex).attitudeIndex = attitudeIndex;
            operatingCases(studyCaseIndex).locationIndex = locationIndex;
            operatingCases(studyCaseIndex).bodyPositionMission = studyBodyPosition;
            operatingCases(studyCaseIndex).R_bodyToMission = studyR;
            operatingCases(studyCaseIndex).cameraPositionBody = studyCameraPositionBody;

            operatingCases(studyCaseIndex).visible = all(cornersCamera(3,:) > 0) && all(insideImage);

            % Plot the projected square in this attitude/altitude panel.
            closedOrder = [1:4, 1];

            plot(ax, studyPixels(1,closedOrder), studyPixels(2,closedOrder), '-o', 'LineWidth', 1);
        end

        axis(ax, 'equal');
        xlim(ax, [0.5, imageWidth + 0.5]);
        ylim(ax, [0.5, imageHeight + 0.5]);
        set(ax, 'YDir', 'reverse');
        grid(ax, 'on');

        xlabel(ax, 'u [px]');
        ylabel(ax, 'v [px]');

        title(ax, sprintf('%s, %.1f m', studyAttitudeNames{attitudeIndex}, altitude_m));
    end
end

sgtitle('Operating-condition geometry check');

studyLegend = legend(ax, studyLocationNames, 'Orientation', 'horizontal');
studyLegend.Layout.Tile = 'south';

%% Report geometric visibility
geometryCoverage = table({operatingCases.name}', [operatingCases.visible]', 'VariableNames', {'Condition', 'AllCornersVisible'});

disp(geometryCoverage);

fprintf('Generated %d operating conditions; %d have all corners visible.\n', numel(operatingCases), sum([operatingCases.visible]));

%% Monte Carlo study: combined camera uncertainty
rng(1);

pixelStdValues = 1; % Corner-coordinate standard deviation [px]
nTrials = 1000;
% Select geometrically valid operating conditions
validCaseIndices = find([operatingCases.visible]);
studyCases = operatingCases(validCaseIndices);

nCases = numel(studyCases);

assert(nCases > 0, 'No operating conditions have fully visible markers.');

fprintf('\nRunning Monte Carlo for %d of %d operating conditions.\n', nCases, numel(operatingCases));
nNoiseLevels = numel(pixelStdValues);

positionReference = 0.10; % Proposed body-frame 3D accuracy reference [m]

%% Provisional calibration uncertainty settings
focalFractionStd = 0.01;       % Fractional focal-length error: 1%
principalPointStd = 10;       % Standard deviation per coordinate [px]
residualDistortionEdgeStd = 5; % Radial displacement std at image corner [px]

% Convert the image-corner displacement allowance into a dimensionless
% first-order radial-distortion coefficient standard deviation.
cornerRadiusPixels = hypot(imageWidth/2, imageHeight/2);

cornerRadiusSquared = (imageWidth/(2*fx))^2 + (imageHeight/(2*fy))^2;

residualK1Std = residualDistortionEdgeStd / (cornerRadiusSquared * cornerRadiusPixels);

%% Generate uncertainty samples
% Corner errors vary by corner and image coordinate.
standardCornerNoise = randn(2, 4, nTrials);

% Calibration errors represent persistent errors in each trial's camera.
% Reuse each trial's calibration draw across the geometry cases.
%
% A common fractional focal error is applied to fx and fy.
% Other parameter groups are sampled independently as an assumption.
focalFractionError_mc = focalFractionStd * randn(1, nTrials);
principalPointError_mc = principalPointStd * randn(2, nTrials);
residualK1_mc = residualK1Std * randn(1, nTrials);

% Dimensions: trial, pixel-noise level, geometry case
bodyErrors_mc = zeros(nTrials, nNoiseLevels, nCases);
missionErrors_mc = zeros(nTrials, nNoiseLevels, nCases);

%% Provisional mounting uncertainty scenario
% Exploratory assumptions, not measured installation accuracy.
mountPositionStd_m = [0.005; 0.005; 0.005]; % Per body axis [m]
mountAngleStd_deg  = [0.5; 0.5; 0.5];       % About body x, y, z [deg]

% Each trial represents one possible mounting-estimation error.
% Reuse these errors across geometry cases for consistent comparisons.
mountPositionError_mc = mountPositionStd_m .* randn(3, nTrials);
mountAngleError_mc = mountAngleStd_deg .* randn(3, nTrials);

%% Load localization results
% Locate the export in the same folder as this marker script.
markerScriptFolder = fileparts(mfilename('fullpath'));

if isempty(markerScriptFolder)
    markerScriptFolder = pwd;
end

navigationData = load(fullfile(markerScriptFolder, 'markerNav_KalmanNormal.mat'), 'markerNavExport');

nav = navigationData.markerNavExport;

assert(strcmp(nav.positionFrame, 'ENU'), 'Navigation position errors must be expressed in ENU.');

%% Select representative navigation samples
% Initially use the entire exported normal-operation time interval.
% These limits can later exclude initialization if necessary.
navigationStartTime_s = nav.time_s(1);
navigationEndTime_s = nav.time_s(end);

eligibleTimeIndices = find(nav.time_s >= navigationStartTime_s & nav.time_s <= navigationEndTime_s);

assert(~isempty(eligibleTimeIndices), 'No navigation samples are inside the selected time interval.');

nNavigationTrials = size(nav.positionErrorMission_m, 3);

% Sample localization trials and times with replacement.
% Each marker trial receives one matched position/attitude error pair.
navigationTrialIndex_mc = randi(nNavigationTrials, 1, nTrials);

navigationTimeIndex_mc = eligibleTimeIndices(randi(numel(eligibleTimeIndices), 1, nTrials));

aircraftPositionErrorMission_mc = zeros(3, nTrials);
attitudeErrorAngles_mc = zeros(3, nTrials);

for trial = 1:nTrials

    navigationTrial = navigationTrialIndex_mc(trial);
    navigationTimeIndex = navigationTimeIndex_mc(trial);

    aircraftPositionErrorMission_mc(:,trial) = nav.positionErrorMission_m(:, navigationTimeIndex, navigationTrial);

    % The export uses radians; our attitudeErrorRotation helper uses degrees.
    attitudeErrorAngles_mc(:,trial) = rad2deg(nav.attitudeErrorAngles_rad(:,navigationTrial));
end

%% Marker-size knowledge uncertainty
% Provisional sensitivity assumption, not verified printing accuracy.
markerSideLengthStd_m = 0.002; % Standard deviation of side-length error [m]

% The true marker uses markerSideLength.
% The pose solver uses an imperfectly known side length.
% Each trial's size error is shared across the geometry cases.
markerSideLengthForSolver_mc = markerSideLength + markerSideLengthStd_m * randn(1, nTrials);

assert(all(markerSideLengthForSolver_mc > 0), 'The assumed marker side lengths must be positive.');

%% Run combined uncertainty simulation
for caseIndex = 1:nCases

    % True aircraft pose and mounting for this operating condition.
    bodyPositionMission = studyCases(caseIndex).bodyPositionMission;
    R_bodyToMission = studyCases(caseIndex).R_bodyToMission;
    cameraPositionBody = studyCases(caseIndex).cameraPositionBody;

    markerCornersBody = R_bodyToMission' * (markerCornersMission - bodyPositionMission);

    markerCornersCamera = R_cameraToBody' * (markerCornersBody - cameraPositionBody);

    assert(all(markerCornersCamera(3,:) > 0), 'Marker corners must be in front of the camera.');

    % Perfect pixel coordinates using the nominal true camera.
    idealPixels = [ fx * markerCornersCamera(1,:) ./ markerCornersCamera(3,:) + cx;
                    fy * markerCornersCamera(2,:) ./ markerCornersCamera(3,:) + cy];

    insideImage = (idealPixels(1,:) >= 0.5 & idealPixels(1,:) <= imageWidth + 0.5) & (idealPixels(2,:) >= 0.5 & idealPixels(2,:) <= imageHeight + 0.5);

    assert(all(insideImage), 'Ideal marker corners are outside the image in case %d.', caseIndex);

    markerCenterBody_true = R_bodyToMission' * (markerCenterMission - bodyPositionMission);

    % Normalized image coordinates: dimensionless ratios Xc/Zc and Yc/Zc.
    xNormalized = (idealPixels(1,:) - cx) / fx;
    yNormalized = (idealPixels(2,:) - cy) / fy;
    radiusSquared = xNormalized.^2 + yNormalized.^2;

    for noiseIndex = 1:nNoiseLevels

        pixelStd = pixelStdValues(noiseIndex);

        for trial = 1:nTrials

            % Model residual radial distortion remaining AFTER correction.
            % One coefficient produces a smooth error field across the image.
            radialScale = 1 + residualK1_mc(trial) * radiusSquared;

            pixelsAfterCorrection = [ fx * (xNormalized .* radialScale) + cx;
                                      fy * (yNormalized .* radialScale) + cy];

            % Add corner-localization error to the corrected pixel coordinates.
            measuredPixels = pixelsAfterCorrection + pixelStd * standardCornerNoise(:,:,trial);

            % Model an imperfect calibration used by the pose estimator.
            estimatedFocalLengths = [fx, fy] * (1 + focalFractionError_mc(trial));

            estimatedPrincipalPoint = [cx, cy] + principalPointError_mc(:,trial)';

            intrinsicsForTrial = cameraIntrinsics(estimatedFocalLengths, estimatedPrincipalPoint, [imageHeight, imageWidth]);

            % Residual distortion is unknown to the estimator.
            % Therefore, no additional distortion correction is applied here.
            %% Construct the marker dimensions assumed by the pose solver
            markerSizeScale = markerSideLengthForSolver_mc(trial) / markerSideLength;

            markerCornersLocal_est = markerCornersLocal * markerSizeScale;

            % Estimate pose using imperfect knowledge of the physical marker size.
            markerToCamera = estimateExtrinsics(measuredPixels', markerCornersLocal_est, intrinsicsForTrial);

            centerCamera_est = markerToCamera.Translation';

            % Stop for investigation if the solver produces an invalid pose.
            estimatedCornersCamera = markerToCamera.R * [markerCornersLocal_est'; zeros(1,4)] + centerCamera_est;

            assert(all(isfinite(estimatedCornersCamera(:))) && all(estimatedCornersCamera(3,:) > 0), ...
                'Invalid pose: case %d, noise level %d, trial %d.', caseIndex, noiseIndex, trial);

            %% Construct this trial's imperfect mounting estimate

            % Estimated camera-origin position relative to the body origin [m].
            cameraPositionBody_est = cameraPositionBody + mountPositionError_mc(:,trial);

            % Mounting-orientation estimation errors about the body axes [deg].
            angleErrorX = mountAngleError_mc(1,trial);
            angleErrorY = mountAngleError_mc(2,trial);
            angleErrorZ = mountAngleError_mc(3,trial);

            Rx_error = [1, 0, 0;
                        0, cosd(angleErrorX), -sind(angleErrorX);
                        0, sind(angleErrorX),  cosd(angleErrorX)];

            Ry_error = [ cosd(angleErrorY), 0, sind(angleErrorY);
                        0, 1, 0;
                        -sind(angleErrorY), 0, cosd(angleErrorY)];

            Rz_error = [cosd(angleErrorZ), -sind(angleErrorZ), 0;
                        sind(angleErrorZ),  cosd(angleErrorZ), 0;
                        0, 0, 1];

            % Apply the alignment error in body coordinates.
            R_mountError = Rz_error * Ry_error * Rx_error;
            R_cameraToBody_est = R_mountError * R_cameraToBody;

            %% Recover coordinates using the imperfect mounting estimate
            centerBody_est = cameraPositionBody_est + R_cameraToBody_est * centerCamera_est;

            %% Aircraft pose estimate from the sampled localization errors

            % Estimated aircraft position in mission ENU coordinates [m].
            bodyPositionMission_est = bodyPositionMission + aircraftPositionErrorMission_mc(:,trial);

            % Estimated attitude using the matching localization-trial error.
            R_attitudeError = attitudeErrorRotation(attitudeErrorAngles_mc(:,trial));

            R_bodyToMission_est = R_bodyToMission * R_attitudeError;

            % Recover marker coordinates using estimated aircraft position and attitude.
            centerMission_est = bodyPositionMission_est + R_bodyToMission_est * centerBody_est;

            % Store 3D error magnitudes [m].
            bodyErrors_mc(trial,noiseIndex,caseIndex) = norm(centerBody_est - markerCenterBody_true);

            missionErrors_mc(trial,noiseIndex,caseIndex) = norm(centerMission_est - markerCenterMission);
        end
    end
end

%% Summarize combined operating-condition results
% This summary assumes one corner-noise setting, as currently configured.
assert(size(bodyErrors_mc, 2) == 1, 'This summary expects one corner-noise setting.');

nSamples = size(bodyErrors_mc, 1);

% Rows = Monte Carlo trials; columns = operating conditions.
bodyErrors = reshape(bodyErrors_mc(:,1,:), nSamples, nCases);
missionErrors = reshape(missionErrors_mc(:,1,:), nSamples, nCases);

bodyRMSE_cm = 100 * sqrt(mean(bodyErrors.^2, 1))';
missionRMSE_cm = 100 * sqrt(mean(missionErrors.^2, 1))';

sortedBodyErrors = sort(bodyErrors, 1);
sortedMissionErrors = sort(missionErrors, 1);

percentileIndex = ceil(0.95 * nSamples);

bodyP95_cm = 100 * sortedBodyErrors(percentileIndex,:)';
missionP95_cm = 100 * sortedMissionErrors(percentileIndex,:)';

bodyWithin_pct = 100 * mean(bodyErrors <= positionReference, 1)';
missionWithin_pct = 100 * mean(missionErrors <= positionReference, 1)';

%% Build a table with one row per operating condition
caseAltitude_m = [studyCases.altitude_m]';

caseAttitude = string(studyAttitudeNames([studyCases.attitudeIndex]))';
caseAttitude = caseAttitude(:);

locationIndices = [studyCases.locationIndex]';

caseImageLocation = compose('%.2f, %.2f', studyImageFractions(locationIndices,1), studyImageFractions(locationIndices,2));

operatingResults = table( ...
    (1:nCases)', ...
    caseAltitude_m, ...
    caseAttitude, ...
    caseImageLocation, ...
    bodyRMSE_cm, ...
    bodyP95_cm, ...
    bodyWithin_pct, ...
    missionRMSE_cm, ...
    missionP95_cm, ...
    missionWithin_pct, ...
    'VariableNames', { ...
    'CaseID', 'Altitude_m', 'Attitude', 'ImageFractions_uv', ...
    'BodyRMSE_cm', 'BodyP95_cm', 'BodyWithinReference_pct', ...
    'MissionRMSE_cm', 'MissionP95_cm', 'MissionWithinReference_pct'});

% Keep the full table in the workspace.
% Print the ten conditions with the largest mission-frame P95.
rankedResults = sortrows(operatingResults, 'MissionP95_cm', 'descend');

fprintf('\nConditions with the largest mission-frame P95:\n');
disp(rankedResults(1:min(10,height(rankedResults)),:));

%% Arrange P95 values into grids for plotting

% Dimensions: altitude, image location, attitude, coordinate frame.
% Frame 1 = body; frame 2 = mission.
% NaN leaves unavailable conditions visibly blank.
p95Grid = nan(nAltitudes, nLocations, nAttitudes, 2);

for caseIndex = 1:nCases

    altitudeIndex = find(studyAltitudes_m == studyCases(caseIndex).altitude_m, 1);

    attitudeIndex = studyCases(caseIndex).attitudeIndex;
    locationIndex = studyCases(caseIndex).locationIndex;

    p95Grid(altitudeIndex,locationIndex,attitudeIndex,1) = bodyP95_cm(caseIndex);

    p95Grid(altitudeIndex,locationIndex,attitudeIndex,2) = missionP95_cm(caseIndex);
end

%% Plot body- and mission-frame P95 heatmaps
locationLabels = compose('%.2f, %.2f', studyImageFractions(:,1), studyImageFractions(:,2));

frameNames = {'Body frame', 'Mission frame'};

% Use one color scale across every panel.
colorMaximum = max([bodyP95_cm; missionP95_cm; 100*positionReference]);

tiledlayout(nAttitudes, 2,'TileSpacing', 'compact', 'Padding', 'compact');

for attitudeIndex = 1:nAttitudes
    for frameIndex = 1:2

        ax = nexttile;
        values = p95Grid(:,:,attitudeIndex,frameIndex);

        imagesc(ax, values, 'AlphaData', isfinite(values));
        hold(ax, 'on');

        set(ax, ...
            'YDir', 'normal', ...
            'Color', [0.85, 0.85, 0.85], ...
            'XTick', 1:nLocations, ...
            'XTickLabel', locationLabels, ...
            'YTick', 1:nAltitudes, ...
            'YTickLabel', compose('%.2g', studyAltitudes_m(:)));

        xlim(ax, [0.5, nLocations + 0.5]);
        ylim(ax, [0.5, nAltitudes + 0.5]);
        xtickangle(ax, 45);

        colormap(ax, parula(256));
        clim(ax, [0, colorMaximum]);

        cb = colorbar(ax);
        cb.Label.String = 'P95 3D error [cm]';

        % Label each available cell and highlight reference exceedances.
        for altitudeIndex = 1:nAltitudes
            for locationIndex = 1:nLocations

                value = values(altitudeIndex,locationIndex);

                if isfinite(value)
                    text(ax, locationIndex, altitudeIndex, ...
                        sprintf('%.1f', value), ...
                        'HorizontalAlignment', 'center', ...
                        'FontSize', 8, ...
                        'Color', 'k', ...
                        'BackgroundColor', 'w', ...
                        'Margin', 0.5);

                    if value > 100*positionReference
                        rectangle(ax, 'Position', ...
                            [locationIndex-0.5, altitudeIndex-0.5, 1, 1], ...
                            'EdgeColor', 'r', 'LineWidth', 1.8);
                    end
                end
            end
        end

        xlabel(ax, 'Marker-center image fractions (u, v)');
        ylabel(ax, 'Body-origin altitude [m]');

        title(ax, sprintf('%s — %s', studyAttitudeNames{attitudeIndex}, frameNames{frameIndex}));
    end
end

sgtitle(sprintf('3D P95 error — red borders exceed the %.0f cm comparison reference', 100*positionReference));

%% Functions
function R_bodyToMission = bodyToMissionFromAngles(rollDeg, pitchDeg, yawDeg)

    Rx = [1, 0, 0;
          0, cosd(rollDeg), -sind(rollDeg);
          0, sind(rollDeg),  cosd(rollDeg)];

    Ry = [ cosd(pitchDeg), 0, sind(pitchDeg);
           0,             1, 0;
          -sind(pitchDeg), 0, cosd(pitchDeg)];

    Rz = [cosd(yawDeg), -sind(yawDeg), 0;
          sind(yawDeg),  cosd(yawDeg), 0;
          0,             0,            1];

    % Level, east-facing FRD body axes expressed in mission ENU.
    R_levelEast = diag([1, -1, -1]);

    R_bodyToMission = R_levelEast * Rz * Ry * Rx;
end

function R_error = attitudeErrorRotation(errorAnglesDeg)

    ex = errorAnglesDeg(1);
    ey = errorAnglesDeg(2);
    ez = errorAnglesDeg(3);

    Rx = [1, 0, 0;
          0, cosd(ex), -sind(ex);
          0, sind(ex),  cosd(ex)];

    Ry = [ cosd(ey), 0, sind(ey);
           0, 1, 0;
          -sind(ey), 0, cosd(ey)];

    Rz = [cosd(ez), -sind(ez), 0;
          sind(ez),  cosd(ez), 0;
          0, 0, 1];

    R_error = Rz * Ry * Rx;
end