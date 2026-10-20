function [p_direct, v_direct] = runDirectRtkImu(a_imu_mission, p_rtk, v_rtk, rtkIndices, dt, rtkAvailable)

    N = size(a_imu_mission, 2);

    p_direct = zeros(3, N);
    v_direct = zeros(3, N);

    % Initialize using the first RTK measurement at t = 0
    p_direct(:, 1) = p_rtk(:, 1);
    v_direct(:, 1) = v_rtk(:, 1);

    nextRtk = 2;

    for k = 2:N

        % Acceleration used over the previous time interval
        acceleration = a_imu_mission(:, k-1);

        % Predict position and velocity using IMU acceleration
        p_direct(:, k) = p_direct(:, k-1) + v_direct(:, k-1) * dt + 0.5 * acceleration * dt^2;

        v_direct(:, k) = v_direct(:, k-1) + acceleration * dt;

        % Check whether this is an RTK measurement time
        if nextRtk <= length(rtkIndices) && k == rtkIndices(nextRtk)

            % Reset position and velocity only if the measurement is available
            if rtkAvailable(nextRtk)

                p_direct(:, k) = p_rtk(:, nextRtk);
                v_direct(:, k) = v_rtk(:, nextRtk);

            end

            % Advance past this measurement, including unavailable measurements
            nextRtk = nextRtk + 1;

        end
    end
end

