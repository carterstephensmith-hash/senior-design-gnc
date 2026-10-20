function [t, p_true, v_true, a_true] = sampleLawnmowerPath(waypoints, maxSpeed, maxAcceleration, dt)
    % Calculate the duration of each segment
    nSegments = size(waypoints,1) - 1;
    segmentDurations = zeros(1,nSegments);

    for i = 1:nSegments
        segmentVector = waypoints(i+1,:) - waypoints(i,:);
        segmentLength = norm(segmentVector);

        [~, ~, segmentDurations(i), ~] = calculateSegmentTiming(segmentLength, maxSpeed, maxAcceleration);
    end

    % Calculate segment start and finish times
    segmentEndTimes = cumsum(segmentDurations);
    segmentStartTimes = [0, segmentEndTimes(1:end-1)];

    totalDuration = segmentEndTimes(end);

    % Use equally spaced samples; the final sample is at or after arrival
    nSteps = ceil(totalDuration/dt);
    t = (0:nSteps)*dt;

    % Allocate the 3D motion histories
    N = length(t);
    p_true = zeros(3,N);
    v_true = zeros(3,N);
    a_true = zeros(3,N);

    % Evaluate the motion at every simulation time
    segmentIndex = 1;

    for k = 1:N

        % Advance to the next segment when necessary
        while segmentIndex < nSegments && t(k) >= segmentEndTimes(segmentIndex)

            segmentIndex = segmentIndex + 1;
        end

        localTime = t(k) - segmentStartTimes(segmentIndex);

        [p_true(:,k), v_true(:,k), a_true(:,k)] = getSegmentState(waypoints(segmentIndex,:), waypoints(segmentIndex+1,:), localTime, maxSpeed, maxAcceleration);
    end
end

%% Function to get accurate path timing accounting for accleration and deaccleration
function [tAccel, tCruise, tTotal, peakSpeed] = calculateSegmentTiming(segmentLength, maxSpeed, maxAcceleration)
    % Time and distance needed to reach maximum speed
    timeToMaxSpeed = maxSpeed/maxAcceleration;
    distanceToMaxSpeed = 0.5*maxAcceleration*timeToMaxSpeed^2;

    % Allow the same distance for braking
    if segmentLength >= 2*distanceToMaxSpeed

        % Long segment: accelerate, cruise, then brake
        peakSpeed = maxSpeed;
        tAccel = timeToMaxSpeed;

        cruiseDistance = segmentLength - 2*distanceToMaxSpeed;
        tCruise = cruiseDistance/peakSpeed;

    else

        % Short segment: accelerate, then immediately brake
        tAccel = sqrt(segmentLength/maxAcceleration);
        peakSpeed = maxAcceleration*tAccel;
        tCruise = 0;

    end
    tTotal = 2*tAccel + tCruise;
end

function [position, velocity, acceleration] = getSegmentState(startPoint, endPoint, localTime, maxSpeed, maxAcceleration)
    % Ensure the endpoints are column vectors
    startPoint = startPoint(:);
    endPoint = endPoint(:);
    
    % Segment geometry
    segmentVector = endPoint - startPoint;
    segmentLength = norm(segmentVector);
    direction = segmentVector/segmentLength;

    % Segment timing
    [tAccel, tCruise, tTotal, peakSpeed] = calculateSegmentTiming(segmentLength, maxSpeed, maxAcceleration);

    % Motion along the segment
    if localTime < tAccel

        distance = 0.5*maxAcceleration*localTime^2;
        speed = maxAcceleration*localTime;
        alongAccel = maxAcceleration;

    elseif localTime < tAccel + tCruise

        accelDistance = 0.5*maxAcceleration*tAccel^2;
        cruiseTime = localTime - tAccel;

        distance = accelDistance + peakSpeed*cruiseTime;
        speed = peakSpeed;
        alongAccel = 0;

    elseif localTime < tTotal

        remainingTime = tTotal - localTime;

        distance = segmentLength - 0.5*maxAcceleration*remainingTime^2;
        speed = maxAcceleration*remainingTime;
        alongAccel = -maxAcceleration;

    else

        distance = segmentLength;
        speed = 0;
        alongAccel = 0;

    end
    % Convert motion along the segment into 3D vectors
    position = startPoint + direction*distance;
    velocity = direction*speed;
    acceleration = direction*alongAccel;
end

