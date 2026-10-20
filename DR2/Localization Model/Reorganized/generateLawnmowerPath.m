%% Waypoint Function
function waypoints = generateLawnmowerPath(areaWidth, areaLength, trackSpacing, altitude)

    % Locate the parallel search passes
    xLines = 0:trackSpacing:areaWidth;

    % Include a pass at the far edge
    if xLines(end) < areaWidth
        xLines(end+1) = areaWidth;
    end

    % Allocate two [x, y, z] endpoints per pass
    waypoints = zeros(2*length(xLines), 3);

    for i = 1:length(xLines)
        x = xLines(i);

        % Alternate the direction of neighboring passes
        if mod(i,2) == 1
            yStart = 0;
            yEnd = areaLength;
        else
            yStart = areaLength;
            yEnd = 0;
        end

        waypoints(2*i-1,:) = [x, yStart, altitude];
        waypoints(2*i,:)   = [x, yEnd, altitude];
    end
end

