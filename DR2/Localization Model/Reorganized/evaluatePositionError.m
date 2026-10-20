function [positionError, positionRMSE, maxPositionError] = evaluatePositionError(p_est, p_true)

    % Position error components: East, North, Up
    positionDifference = p_est - p_true;

    % Total 3D position error at each time step
    positionError = sqrt(sum(positionDifference.^2, 1));

    % Summarize error over the flight
    positionRMSE = sqrt(mean(positionError.^2));
    maxPositionError = max(positionError);

end

