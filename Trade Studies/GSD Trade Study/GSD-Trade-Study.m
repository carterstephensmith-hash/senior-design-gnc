clear; clc; close all;

%% GSD Trade Study Inputs
% Generic image-size candidates, not specific camera recommendations.
hFOVValues = [60, 80, 100]; % degrees
imageResolutions = [1280, 720;  % 720p
                    1920, 1080; % 1080p
                    2560, 1440; % 2k
                    3840, 2160]; % 4k  % [width, height] in pixels

markerSize = 0.08;          % m, square marker side length
minimumMarkerPixels = 60;   % required pixels on each side
maximumAltitude = 2;       % m, mission altitude limit

%% Run Study
maximumGSD = markerSize / minimumMarkerPixels;
fprintf('Maximum allowable GSD: %.4f cm/pixel\n', 100 * maximumGSD);

gsdResults = runGSDTrade(hFOVValues, imageResolutions, markerSize, minimumMarkerPixels, maximumAltitude);

disp(gsdResults);

%% Local Function
function results = runGSDTrade(hFOVValues, resolutions, markerSize, minimumMarkerPixels, maximumAltitude)

    maximumGSD = markerSize / minimumMarkerPixels; % m/pixel
    numberOfResolutions = size(resolutions, 1);
    numberOfFOVs = numel(hFOVValues);

    data = zeros(numberOfResolutions * numberOfFOVs, 10);
    allowableAltitudes = zeros(numberOfResolutions, numberOfFOVs);
    row = 0;

    for i = 1:numberOfResolutions
        imageWidth = resolutions(i, 1);
        imageHeight = resolutions(i, 2);

        for j = 1:numberOfFOVs
            hFOV = hFOVValues(j);
            vFOV = 2 * atand((imageHeight / imageWidth) * tand(hFOV / 2));

            % Solve for the highest altitude meeting both pixel limits.
            horizontalAltitudeLimit = maximumGSD * imageWidth / (2 * tand(hFOV / 2));
            verticalAltitudeLimit = maximumGSD * imageHeight / (2 * tand(vFOV / 2));

            gsdAltitudeLimit = min(horizontalAltitudeLimit, verticalAltitudeLimit);
            altitude = min(gsdAltitudeLimit, maximumAltitude);

            % Evaluate the footprint and marker at this usable altitude.
            footprintWidth = 2 * altitude * tand(hFOV / 2);
            footprintLength = 2 * altitude * tand(vFOV / 2);
            hGSD = footprintWidth / imageWidth;
            vGSD = footprintLength / imageHeight;
            horizontalPixels = markerSize / hGSD;
            verticalPixels = markerSize / vGSD;

            allowableAltitudes(i, j) = altitude;
            row = row + 1;
            data(row, :) = [imageWidth, imageHeight, hFOV, gsdAltitudeLimit, altitude, footprintWidth, footprintLength, ...
                100 * max(hGSD, vGSD), horizontalPixels, verticalPixels];
        end
    end

    results = array2table(data, 'VariableNames', {'ImageWidth_px', 'ImageHeight_px', 'HorizontalFOV_deg', 'GSDAltitudeLimit_m', 'MaximumUsableAltitude_m', ...
        'FootprintWidth_m', 'FootprintLength_m', 'GSD_cm_perPixel', 'HorizontalMarkerPixels', 'VerticalMarkerPixels'});

    %% Plot Maximum Usable Altitude
    figure('Color', 'w', 'WindowStyle', 'normal', 'Name', 'GSD Trade Study');
    hold on;
    colors = lines(numberOfResolutions);

    for i = 1:numberOfResolutions
        plot(hFOVValues, allowableAltitudes(i, :), '-o', 'Color', colors(i, :), 'LineWidth', 1.5, 'MarkerSize', 7, ...
            'DisplayName', sprintf('%d x %d', resolutions(i, 1), resolutions(i, 2)));
    end

    yline(maximumAltitude, '--', 'Mission altitude limit', 'HandleVisibility', 'off');
    xlabel('Horizontal FOV (degrees)');
    ylabel('Maximum usable altitude (m)');
    title({sprintf('GSD Trade: At Least %d x %d Marker Pixels', minimumMarkerPixels, minimumMarkerPixels), 'At each tested FOV, altitudes at or below its point are feasible'});
    xticks(hFOVValues);
    ylim([0, maximumAltitude * 1.15]);
    legend('Location', 'best');
    grid on;
    box on;
end
