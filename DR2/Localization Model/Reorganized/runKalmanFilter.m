function [p_kf, v_kf, b_kf] = runKalmanFilter(a_imu_mission, f_imu_body, R_est_bodyToMission, p_rtk, v_rtk, rtkIndices, dt, ...
    x0_kf, P0_kf, H, R_rtk, accelNoiseStd, accelBiasStd, biasTau, rtkAvailable)

    N = size(a_imu_mission, 2);

    % Store the nine estimated states throughout the flight
    x_history = zeros(12, N);
    x_history(:, 1) = x0_kf;

    % Current state estimate and covariance
    x = x0_kf;
    P = P0_kf;

    nextRtk = 2;

    % This contribution is constant for our chosen dt and bias model
    Q_bias = buildBiasProcessNoise(dt, accelBiasStd, biasTau);

    for k = 2:N

        % Inputs for the interval from k-1 to k
        acceleration = a_imu_mission(:, k-1);
        specificForceBody = f_imu_body(:, k-1);
        R_est = R_est_bodyToMission(:, :, k-1);

        % Predict the state
        x_pred = predictKalmanState(x, acceleration, specificForceBody, R_est, dt, biasTau);

        % Predict the covariance
        F = buildStateTransition(R_est, specificForceBody, dt, biasTau);
        Q_accel = buildAccelProcessNoise(R_est, dt, accelNoiseStd);
        Q = Q_accel + Q_bias;

        P_pred = F * P * F' + Q;

        % Retain the prediction unless an available RTK measurement corrects it
        x = x_pred;
        P = P_pred;

        if nextRtk <= length(rtkIndices)  && k == rtkIndices(nextRtk)

            if rtkAvailable(nextRtk)

                z_rtk = [p_rtk(:, nextRtk);
                        v_rtk(:, nextRtk)];

                [x, P] = updateWithRtk(x_pred, P_pred, z_rtk, H, R_rtk);

            end

            % Advance even when this measurement is unavailable
            nextRtk = nextRtk + 1;

        end

        x_history(:, k) = x;

    end

    % Separate the state history into useful outputs
    p_kf = x_history(1:3, :);
    v_kf = x_history(4:6, :);
    b_kf = x_history(7:9, :);

end

function x_pred = predictKalmanState(x_previous, accelMissionSample, specificForceBody, R_est, dt, biasTau)

    % Extract the previous estimates
    position = x_previous(1:3);
    velocity = x_previous(4:6);
    biasBody = x_previous(7:9);
    attitudeErrorEstimate = x_previous(10:12);

    % Express the estimated bias in mission coordinates
    biasMission = R_est * biasBody;

    % Acceleration correction for estimated attitude error
    attitudeSensitivity = buildAttitudeSensitivity(R_est, specificForceBody);

    attitudeCorrection =  attitudeSensitivity * attitudeErrorEstimate;

    % Correct the measured acceleration using the estimated bias
    correctedAcceleration = accelMissionSample - biasMission + attitudeCorrection;

    % Predict position and velocity
    positionPred = position + velocity * dt + 0.5 * correctedAcceleration * dt^2;

    velocityPred = velocity + correctedAcceleration * dt;

    % Predict bias using our Gauss–Markov model
    alpha = exp(-dt / biasTau);
    biasPred = alpha * biasBody;

    % Our modeled attitude offset persists throughout the flight
    attitudeErrorPred = attitudeErrorEstimate;

    % Reassemble the nine-state vector
    x_pred = [positionPred;
              velocityPred;
              biasPred;
              attitudeErrorPred];

end

function F = buildStateTransition(R_est, specificForceBody, dt, biasTau)

    I3 = eye(3);
    Z3 = zeros(3);

    alpha = exp(-dt / biasTau);

    attitudeSensitivity = buildAttitudeSensitivity(R_est, specificForceBody);

    F = [I3, dt*I3, -0.5*dt^2*R_est, 0.5*dt^2*attitudeSensitivity;
         Z3,    I3,       -dt*R_est,       dt*attitudeSensitivity;
         Z3,    Z3,       alpha*I3,                          Z3;
         Z3,    Z3,             Z3,                          I3];

end

function Q_accel = buildAccelProcessNoise(R_est, dt, accelNoiseStd)

    % Accelerometer noise covariance in body axes
    accelNoiseCovBody = diag(accelNoiseStd.^2);

    % Map acceleration noise into the nine-state prediction
    G = [0.5 * dt^2 * R_est;
         dt * R_est;
         zeros(3);
         zeros(3)];

    % Contribution to prediction uncertainty
    Q_accel = G * accelNoiseCovBody * G';

end

function Q_bias = buildBiasProcessNoise(dt, accelBiasStd, biasTau)

    alpha = exp(-dt / biasTau);

    % Same random-step standard deviation used in the sensor model
    biasStepStd = accelBiasStd * sqrt(1 - alpha^2);

    Q_bias = zeros(12, 12);

    % New uncertainty directly enters the three bias states
    Q_bias(7:9, 7:9) = biasStepStd^2 * eye(3);

end

function [x_updated, P_updated] = updateWithRtk(x_pred, P_pred, z_rtk, H, R_rtk)

    % Difference between the measurement and predicted measurement
    innovation = z_rtk - H * x_pred;

    % Innovation covariance
    S = H * P_pred * H' + R_rtk;

    % Kalman gain
    K = (P_pred * H') / S;

    % Correct position, velocity, and bias estimates
    x_updated = x_pred + K * innovation;

    % Update the uncertainty of the corrected state
    A = eye(12) - K * H;
    P_updated = A * P_pred * A' + K * R_rtk * K';
end

function attitudeSensitivity = buildAttitudeSensitivity(R_est, specificForceBody)

    f = specificForceBody;

    % Matrix representing cross(f, angleError)
    forceCrossMatrix = [    0, -f(3),  f(2);
                        f(3),     0, -f(1);
                       -f(2),  f(1),     0];

    % Express the acceleration correction in mission coordinates
    attitudeSensitivity = R_est * forceCrossMatrix;

end

