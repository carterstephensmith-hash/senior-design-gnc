clear; clc; close all;

%% Baseline Conditions
baseline.areaWidth = 5;
baseline.areaLength = 5;
baseline.speed = 0.5;

baseline.hFOV = 60;
baseline.altitude = 2;
baseline.overlap = 0.20;

baseline.imageResolution = [1280, 720];
baseline.markerSize = 0.08;
baseline.minimumMarkerPixels = 10;

baseline.sampleSpacing = 0.05;
baseline.gridSpacing = 0.05;
baseline.markerGridSpacing = 0.10;

algorithms = [
    "Lawnmower", "Contracting Square", "Diagonal Lawnmower"
];

fovResults = runOneFactorSweep(algorithms, baseline, "hFOV", [60, 80, 100]);

altitudeResults = runOneFactorSweep(algorithms, baseline, "altitude", [0.5, 1.0, 1.5, 2.0]);

overlapResults = runOneFactorSweep(algorithms, baseline, "overlap", [0.10, 0.20, 0.30, 0.40]);

%% Plot Baseline Algorithm Paths
plotCoveragePaths(algorithms, baseline);

plotSweepResults(fovResults, algorithms);
plotSweepResults(altitudeResults, algorithms);
plotSweepResults(overlapResults, algorithms);

%% Functions

function waypoints = generateLawnmowerPath(area_width, area_length, trackSpacing, altitude)

    xlines = 0:trackSpacing:area_width;
    
    % if the area_width is not divisible by area_length
    if mod(area_width(end), trackSpacing) ~= 0
        xlines = [xlines, area_width];
    end

    waypoints = [];
    
    for i = 1:length(xlines)
        x = xlines(i);

        % Alternate the direction of each pass
        if mod(i,2) == 1
            yStart = 0;
            yEnd = area_length;
        else
            yStart = area_length;
            yEnd = 0;
        end
    
        waypoints = [waypoints;
                     x, yStart, altitude;
                     x, yEnd, altitude];
    
    end

end

function waypoints = generateContractingSquarePath(area_width, area_length, xTrackSpacing, yTrackSpacing, altitude)

    centerX = area_width/2;
    centerY = area_length/2;

    % Initial outer boundaries
    left = 0;
    right = area_width;
    bottom = 0;
    top = area_length;

    % Begin at the lower-left corner
    waypoints = [left, bottom, altitude];

    while true

        % Follow the current contracting rectangle
        waypoints = [waypoints;
                     left,  top,    altitude;
                     right, top,    altitude;
                     right, bottom, altitude];

        % Shrink each boundary without crossing the center
        nextLeft = min(left + xTrackSpacing, centerX);
        nextRight = max(right - xTrackSpacing, centerX);

        nextBottom = min(bottom + yTrackSpacing, centerY);
        nextTop = max(top - yTrackSpacing, centerY);

        % Move inward toward the next layer
        waypoints(end+1,:) = [nextLeft, bottom, altitude];

        xFinished = nextLeft >= nextRight;
        yFinished = nextBottom >= nextTop;

        if xFinished || yFinished

            % Width finished first:
            % fly through the remaining vertical strip
            if xFinished && ~yFinished
                waypoints(end+1,:) = [centerX, nextTop, altitude];
        
            % Height finished first:
            % fly through the remaining horizontal strip
            elseif yFinished && ~xFinished
                waypoints(end+1,:) = [nextLeft, centerY, altitude];
        
                waypoints(end+1,:) = [nextRight, centerY, altitude];
        
            % Both directions finished together
            else
                waypoints(end+1,:) = [centerX, centerY, altitude];
            end
        
            break
        end

        % Continue with the next inner rectangle
        left = nextLeft;
        right = nextRight;
        bottom = nextBottom;
        top = nextTop;
    end
end

function waypoints = generateDiagonalLawnmowerPath(area_width, area_length, xTrackSpacing, yTrackSpacing, altitude)
    
    maxOffsetSpacing = min(xTrackSpacing/2 + yTrackSpacing, xTrackSpacing + yTrackSpacing/2);

    numberOfIntervals = 2*ceil((area_width + area_length)/(2*maxOffsetSpacing));

    offsets = linspace(0, area_width + area_length, numberOfIntervals + 1);

    offsets = unique([offsets, area_width, area_length]);

    offsets = offsets([true, diff(offsets) > 1e-10]);

    waypoints = [0, 0, altitude];

    for i = 2:length(offsets)
    
        offset = offsets(i);

        % First intersection wih the AOI boundary
        xA = max(0, offset - area_length);
        yA = offset - xA;

        % Second intersection with the AOI boundary
        xB = min(area_width, offset);
        yB = offset - xB;

        pointA = [xA, yA, altitude];
        pointB = [xB, yB, altitude];

        % Alternate the direction of successive passes
        if mod(i,2) == 0
            waypoints = [waypoints; pointA; pointB];
        else
            waypoints = [waypoints; pointB; pointA];
        end

        % Remove consecutive duplicate waypoints
        keepWayPoint = [true; any(abs(diff(waypoints,1,1)) > 1e-10, 2)];

        waypoints = waypoints(keepWayPoint,:);
    end
end

function [positions, time] = samplePath(waypoints, sampleSpacing, speed)
    positions = [];

    for i = 1:size(waypoints,1)-1

        startPoint = waypoints(i,:);
        endPoint = waypoints(i+1,:);

        segmentLength = norm(endPoint - startPoint);

        numberOfPoints = max(2, ceil(segmentLength/sampleSpacing) + 1);

        fraction = linspace(0,1,numberOfPoints)';

        segmentPositions = startPoint + fraction.*(endPoint-startPoint);

        % Avoid repeating the connecting point
        if i > 1
            segmentPositions = segmentPositions(2:end,:);
        end

        positions = [positions; segmentPositions];

    end

    positionChanges = diff(positions);

    segmentDistances = sqrt(sum(positionChanges.^2, 2));

    distance = [0; cumsum(segmentDistances)];

    time = distance/speed;

end

% Calculates the ground rectangle visible beneath the drone
function [xCorners, yCorners, width, length] = calculateCameraFootprint(position, hFOV, vFOV)
    x = position(1);
    y = position(2);
    altitude = position(3);

    width = 2*altitude*tand(hFOV/2);
    length = 2*altitude*tand(vFOV/2);

    xCorners = x + [-width/2, width/2, width/2, -width/2];
    yCorners = y + [-length/2, -length/2, length/2, length/2];
    
end

function [coveragePercent, X, Y, covered] = evaluateCoverage(positions, aeraWidth, areaLength, hFOV, vFOV, gridSpacing)
    xGrid = 0:gridSpacing:aeraWidth;
    yGrid = 0:gridSpacing:areaLength;

    [X,Y] = meshgrid(xGrid, yGrid);
    covered = false(size(X));

    for i = 1:size(positions,1)
        [~, ~, width, length] = calculateCameraFootprint(positions(i,:), hFOV, vFOV);

        insideFootprint = abs(X-positions(i,1)) <= width/2 & abs(Y-positions(i,2)) <= length/2;

        covered = covered | insideFootprint;
    end

    coveragePercent = 100*sum(covered(:))/numel(covered);

end

function [pathLength, numberOfTurns] = calculatePathMetrics(waypoints)

    pathSegments = diff(waypoints(:,1:2));

    segmentLengths = sqrt(sum(pathSegments.^2, 2));

    pathLength = sum(segmentLengths);

    % Remove any 0 length segments
    pathSegments = pathSegments(segmentLengths > 0,:);

    headings = atan2(pathSegments(:,2), pathSegments(:,1));

    headingChanges = diff(headings);

    headingChanges = atan2(sin(headingChanges), cos(headingChanges));

    numberOfTurns = sum(abs(headingChanges) > deg2rad(1));
end

function [detectionPercent, averageDetectionTime, worstDetectionTime, markerLocations] = evaluateDetectionMetrics(positions, time, areaWidth, areaLength, markerGridSpacing, ...
                                                                              markerSize, footprintWidth, footprintLength, imageResolution, minimumMarkerPixels)

    % Create Possible marker locations
    xMarker = 0:markerGridSpacing:areaWidth;
    yMarker = 0:markerGridSpacing:areaLength;

    [Xmarker, Ymarker] = meshgrid(xMarker,yMarker);

    markerLocations = [Xmarker(:), Ymarker(:)];

    numberOfMarkers = size(markerLocations,1);
    detectionTimes = NaN(numberOfMarkers,1);

    hGSD = footprintWidth/imageResolution(1);
    vGSD = footprintLength/imageResolution(2);

    horizontalPixels = markerSize/hGSD;
    verticalPixels = markerSize/vGSD;

    markerLargeEnough = horizontalPixels >= minimumMarkerPixels && verticalPixels >= minimumMarkerPixels;

    if markerLargeEnough

        % Test every possible marker location
        for i = 1:numberOfMarkers
    
            markerX = markerLocations(i,1);
            markerY = markerLocations(i,2);
    
            markerVisible = abs(positions(:,1)-markerX) <= footprintWidth/2 & abs(positions(:,2)-markerY) <= footprintLength/2;
    
            firstDetection = find(markerVisible, 1, 'first');
    
            if ~isempty(firstDetection)
                detectionTimes(i) = time(firstDetection);
            end
        end
    end

    detected = ~isnan(detectionTimes);

    detectionPercent = 100*sum(detected)/numberOfMarkers;

    if any(detected)
        averageDetectionTime = mean(detectionTimes(detected));
    else
        averageDetectionTime = NaN;
    end

    if all(detected)
        worstDetectionTime = max(detectionTimes);
    else
        worstDetectionTime = Inf;
    end

end

function results = runOneFactorSweep(algorithms, baseline, variableName, values)

    results = table;

    for i = 1:length(values)

        parameters = baseline;

        % Change only the selected variable
        parameters.(variableName) = values(i);

        for j = 1:length(algorithms)
            
            algorithm = algorithms(j);

            newRow = runCoverageCase(algorithm, parameters);

            newRow.SweepVariable = variableName;
            newRow.SweepValue = values(i);

            results = [results; newRow];
        end
    end
end

function metrics = runCoverageCase(algorithm, parameters)

    vFOV = 2*atand((parameters.imageResolution(2) / parameters.imageResolution(1)) * tand(parameters.hFOV/2));

    footprintWidth = 2*parameters.altitude * tand(parameters.hFOV/2);
    footprintLength = 2*parameters.altitude * tand(vFOV/2);

    trackSpacing = footprintWidth*(1-parameters.overlap);

    % For sprial algorithm
    yTrackSpacing = footprintLength*(1-parameters.overlap);

    % Generate Selected path
    switch algorithm

        case "Lawnmower"
            waypoints = generateLawnmowerPath(parameters.areaWidth, parameters.areaLength, trackSpacing, parameters.altitude);

        case "Contracting Square"
            waypoints = generateContractingSquarePath(parameters.areaWidth, parameters.areaLength, trackSpacing, yTrackSpacing, parameters.altitude);

        case "Diagonal Lawnmower"
            waypoints = generateDiagonalLawnmowerPath(parameters.areaWidth, parameters.areaLength, trackSpacing, yTrackSpacing, parameters.altitude);

        otherwise
            error("Unkown coverage algorithm")
    end

    [positions, time] = samplePath(waypoints, parameters.sampleSpacing, parameters.speed);

    [coveragePercent, ~, ~, ~] = evaluateCoverage(positions, parameters.areaWidth, parameters.areaLength, parameters.hFOV, vFOV, parameters.gridSpacing);

    [pathLength, numberOfTurns] = calculatePathMetrics(waypoints);

    searchTime = time(end);

    [detectionPercent, averageDetectionTime, worstDetectionTime] = evaluateDetectionMetrics(positions, time, parameters.areaWidth, parameters.areaLength, ...
                parameters.markerGridSpacing, parameters.markerSize, footprintWidth, footprintLength, parameters.imageResolution, parameters.minimumMarkerPixels);
    
    hGSD = footprintWidth/parameters.imageResolution(1);
    vGSD = footprintLength/parameters.imageResolution(2);

    hPixels = parameters.markerSize/hGSD;
    vPixels = parameters.markerSize/vGSD;

    metrics = table(algorithm, parameters.hFOV, vFOV, parameters.altitude, parameters.overlap, footprintWidth, footprintLength, ...
    trackSpacing, coveragePercent, detectionPercent, pathLength, searchTime, averageDetectionTime, worstDetectionTime, numberOfTurns, ...
    hPixels, vPixels, 'VariableNames', {'Algorithm', 'HorizontalFOV_deg', 'VerticalFOV_deg', 'Altitude_m', 'Overlap', 'FootprintWidth_m', ...
    'FootprintLength_m', 'TrackSpacing_m', 'CoveragePercent', 'DetectionPercent', 'PathLength_m', 'SearchTime_s', 'AverageDetectionTime_s', ...
    'WorstDetectionTime_s', 'NumberOfTurns', 'HorizontalMarkerPixels', 'VerticalMarkerPixels'});
end

function plotCoveragePaths(algorithms, parameters)

    % Calculate camera footprint
    vFOV = 2*atand((parameters.imageResolution(2) / parameters.imageResolution(1)) * tand(parameters.hFOV/2));

    footprintWidth = 2*parameters.altitude*tand(parameters.hFOV/2);

    footprintLength = 2*parameters.altitude*tand(vFOV/2);

    % Required spacing in each direction
    xTrackSpacing = footprintWidth*(1-parameters.overlap);

    yTrackSpacing = footprintLength*(1-parameters.overlap);

    numberOfAlgorithms = length(algorithms);

    % Create figure
    figure( 'Color', 'white', 'Name', 'Coverage Algorithm Paths', 'Position', [100, 100, 500*numberOfAlgorithms, 475]);

    tiledlayout(1, numberOfAlgorithms, 'TileSpacing', 'compact', 'Padding', 'compact');

    pathColors = lines(numberOfAlgorithms);

    for i = 1:numberOfAlgorithms

        algorithm = algorithms(i);

        % Generate selected path
        switch algorithm

            case "Lawnmower"
                waypoints = generateLawnmowerPath(parameters.areaWidth,parameters.areaLength, xTrackSpacing, parameters.altitude);

            case "Contracting Square"
                waypoints = generateContractingSquarePath(parameters.areaWidth, parameters.areaLength, xTrackSpacing, yTrackSpacing, parameters.altitude);

            case "Diagonal Lawnmower"
                waypoints = generateDiagonalLawnmowerPath(parameters.areaWidth, parameters.areaLength, xTrackSpacing, yTrackSpacing, parameters.altitude);
    
            otherwise
                error("Unknown coverage algorithm: %s", algorithm);
        end

        % Calculate useful path information
        [pathLength, numberOfTurns] = calculatePathMetrics(waypoints);

        searchTime = pathLength/parameters.speed;

        % Create subplot
        nexttile;
        hold on;

        % Draw the search area
        rectangle('Position', [0, 0, parameters.areaWidth, parameters.areaLength], 'FaceColor', [0.96, 0.96, 0.96], 'EdgeColor', [0.15, 0.15, 0.15], 'LineWidth', 1.5);

        % Draw path and waypoints
        plot(waypoints(:,1), waypoints(:,2), '-o', 'Color', pathColors(i,:), 'LineWidth', 2, 'MarkerSize', 4, 'MarkerFaceColor', pathColors(i,:));

        % Mark starting position
        plot(waypoints(1,1), waypoints(1,2), 'o', 'MarkerSize', 10, 'MarkerFaceColor', [0.15, 0.70, 0.25], 'MarkerEdgeColor', 'black');

        % Mark ending position
        plot(waypoints(end,1), waypoints(end,2), 's', 'MarkerSize', 10, 'MarkerFaceColor', [0.85, 0.20, 0.20], 'MarkerEdgeColor', 'black');

        % Start and end labels
        text(waypoints(1,1) + 0.10, waypoints(1,2) + 0.10, 'Start', 'FontWeight', 'bold');

        text(  waypoints(end,1) + 0.10, waypoints(end,2) + 0.10, 'End', 'FontWeight', 'bold');

        % Plot formatting
        xlabel('X Position (m)');
        ylabel('Y Position (m)');

        title({char(algorithm), sprintf('Length: %.1f m | Time: %.1f s | Turns: %d', pathLength, searchTime, numberOfTurns)});

        axis equal;
        grid on;
        box on;

        plotMargin = 0.35;

        xlim([-plotMargin, parameters.areaWidth + plotMargin]);

        ylim([-plotMargin, parameters.areaLength + plotMargin]);

        hold off;
    end

    sgtitle(sprintf( ...
        ['Baseline Coverage Paths: Altitude = %.1f m, ' 'Horizontal FOV = %.0f degrees, Overlap = %.0f%%'], parameters.altitude, ...
        parameters.hFOV, 100*parameters.overlap), 'FontWeight', 'bold');
end

function plotSweepResults(results, algorithms)

    algorithms = string(algorithms);
    sweepVariable = string(results.SweepVariable(1));

    % Configure the horizontal axis
    switch sweepVariable

        case "hFOV"
            xLabel = 'Horizontal FOV (degrees)';
            figureTitle = 'Horizontal FOV Trade Study';
            xScale = 1;

        case "altitude"
            xLabel = 'Altitude (m)';
            figureTitle = 'Altitude Trade Study';
            xScale = 1;

        case "overlap"
            xLabel = 'Desired Overlap (%)';
            figureTitle = 'Overlap Trade Study';
            xScale = 100;
    end

    % Table columns to plot
    metricNames = {
        'SearchTime_s'
        'WorstDetectionTime_s'
        'AverageDetectionTime_s'
        'NumberOfTurns'
        'DetectionPercent'
    };

    panelTitles = {
        'Search Time'
        'Worst-Case Detection Time'
        'Average Detection Time'
        'Number of Turns'
        'Marker Detection'
    };

    yLabels = {
        'Time (s)'
        'Time (s)'
        'Time (s)'
        'Turns'
        'Detection (%)'
    };

    % Consistent algorithm colors and markers
    numberOfAlgorithms = numel(algorithms);
    colors = lines(numberOfAlgorithms);
    markers = {'o', 's', '^', 'd', 'v', '>'};

    xTicks = unique(results.SweepValue)*xScale;

    figure('Color', 'white', 'Name', figureTitle, 'WindowStyle', 'normal');

    tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

    for k = 1:numel(metricNames)

        % Detection plot spans the final two tiles
        if k == 5
            ax = nexttile([1, 2]);
        else
            ax = nexttile;
        end

        hold(ax, 'on');

        lineHandles = gobjects(numberOfAlgorithms,1);
        anyFiniteValues = false;

        for j = 1:numberOfAlgorithms

            % Select and sort this algorithm's results
            selected = string(results.Algorithm) == algorithms(j);

            algorithmResults = sortrows(results(selected,:), 'SweepValue');

            x = algorithmResults.SweepValue*xScale;
            y = algorithmResults.(metricNames{k});

            % Keep coverage as a feasibility check,
            % even though it is no longer plotted
            feasible = algorithmResults.CoveragePercent >= 100 - 1e-8 & algorithmResults.DetectionPercent >= 100 - 1e-8;

            % Compare detection times only for feasible cases
            if k == 2 || k == 3
                y(~feasible) = NaN;
            end

            % Undefined/infinite values appear as gaps
            y(~isfinite(y)) = NaN;

            anyFiniteValues = anyFiniteValues || any(isfinite(y));

            marker = markers{mod(j-1, numel(markers)) + 1};

            lineHandles(j) = plot(ax, x, y, 'Color', colors(j,:), 'LineStyle', '-', 'Marker', marker, 'LineWidth', 1.5, 'MarkerSize', 7);
        end

        title(ax, panelTitles{k});
        xlabel(ax, xLabel);
        ylabel(ax, yLabels{k});

        xticks(ax, xTicks);
        grid(ax, 'on');
        box(ax, 'on');

        set(ax, 'FontSize', 11, 'LineWidth', 1);

        if k == 5

            yline(ax, 100, '--', 'Color', [0.35, 0.35, 0.35], 'HandleVisibility', 'off');

            ylim(ax, [0, 105]);

        elseif anyFiniteValues

            currentLimits = ylim(ax);
            ylim(ax, [0, currentLimits(2)]);

        else

            text(ax, 0.5, 0.5, 'No feasible cases', 'Units', 'normalized', 'HorizontalAlignment', 'center');
        end

        hold(ax, 'off');
    end

    % Shared legend below all panels
    lgd = legend(ax, lineHandles, cellstr(algorithms), 'Orientation', 'horizontal', 'Interpreter', 'none');

    lgd.Layout.Tile = 'south';

    sgtitle({figureTitle, 'Detection-time gaps indicate infeasible cases'}, 'FontWeight', 'bold');
end