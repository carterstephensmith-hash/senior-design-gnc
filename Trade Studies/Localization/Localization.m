%% Housekeeping 

clc, clear, close all;

%% Set up

dt = 0.01;
t = 0:dt:20;

x_true = 2*cos(0.3*t);
y_true = 2*sin(0.3*t);
z_true = 2 + 0.2*sin(0.5*t);

vx_true = gradient(x_true,dt);
vy_true = gradient(y_true,dt);
vz_true = gradient(z_true,dt);

ax_true = gradient(vx_true,dt);
ay_true = gradient(vy_true,dt);
az_true = gradient(vz_true,dt);

gpsNoiseStd = 1; % GPS noise [m]
x_gps = x_true + gpsNoiseStd*randn(size(x_true));
y_gps = y_true + gpsNoiseStd*randn(size(y_true));
z_gps = z_true + gpsNoiseStd*randn(size(z_true));

gpsNoiseVelStd = 0.1; % GPS noise [m/s]
vx_gps = vx_true + gpsNoiseVelStd*randn(size(vx_true));
vy_gps = vy_true + gpsNoiseVelStd*randn(size(vy_true));
vz_gps = vz_true + gpsNoiseVelStd*randn(size(vz_true));

accelNoiseStd = 0.1; % Accelerometer [m/s^2]
ax_IMU = ax_true + accelNoiseStd*randn(size(t));
ay_IMU = ay_true + accelNoiseStd*randn(size(t));
az_IMU = az_true + accelNoiseStd*randn(size(t));

%% GPS ONLY

x_est_GPS = x_gps;
y_est_GPS = y_gps;
z_est_GPS = z_gps;

vx_est_GPS = vx_gps;
vy_est_GPS = vy_gps;
vz_est_GPS = vz_gps;

PositionError_GPS = sqrt((x_est_GPS-x_true).^2 + (y_est_GPS-y_true).^2 + (z_est_GPS-z_true).^2);
VelError_GPS = sqrt((vx_est_GPS-vx_true).^2 + (vy_est_GPS-vy_true).^2 + (vz_est_GPS-vz_true).^2);

%% IMU using Kalman filter (x_k = F*x_(k-1) + B*a_IMU)

% State vectors
x_est_IMU = zeros(size(t));
y_est_IMU = zeros(size(t));
z_est_IMU = zeros(size(t));
vx_est_IMU = zeros(size(t));
vy_est_IMU = zeros(size(t));
vz_est_IMU = zeros(size(t));

% Initial conditons 
x_est_IMU(1) = x_true(1);
y_est_IMU(1) = y_true(1);
z_est_IMU(1) = z_true(1);
vx_est_IMU(1) = vx_true(1);
vy_est_IMU(1) = vy_true(1);
vz_est_IMU(1) = vz_true(1);

% Kalman filter
for k = 2:length(t)
    
    % Position (reduces to kinematics)
    x_est_IMU(k) = x_est_IMU(k-1) + vx_est_IMU(k-1)*dt + 0.5*ax_IMU(k)*dt^2;
    y_est_IMU(k) = y_est_IMU(k-1) + vy_est_IMU(k-1)*dt + 0.5*ay_IMU(k)*dt^2;
    z_est_IMU(k) = z_est_IMU(k-1) + vz_est_IMU(k-1)*dt + 0.5*az_IMU(k)*dt^2;

    % Velocity (reduces to kinematics)
    vx_est_IMU(k) = vx_est_IMU(k-1) + ax_IMU(k)*dt;
    vy_est_IMU(k) = vy_est_IMU(k-1) + ay_IMU(k)*dt;
    vz_est_IMU(k) = vz_est_IMU(k-1) + az_IMU(k)*dt;

end

PositionError_IMU = sqrt((x_est_IMU-x_true).^2 + (y_est_IMU-y_true).^2 + (z_est_IMU-z_true).^2);
VelError_IMU = sqrt((vx_est_IMU-vx_true).^2 + (vy_est_IMU-vy_true).^2 + (vz_est_IMU-vz_true).^2);


%% GPS + IMU

% State vector
% x = [x; y; z; vx; vy; vz]

% State transition matrix
F = [1 0 0 dt 0 0;
     0 1 0 0 dt 0;
     0 0 1 0 0 dt;
     0 0 0 1 0 0;
     0 0 0 0 1 0;
     0 0 0 0 0 1];

% IMU acceleration input matrix
B = [0.5*dt^2 0 0;
     0 0.5*dt^2 0;
     0 0 0.5*dt^2;
     dt 0 0;
     0 dt 0;
     0 0 dt];

% GPS measures position AND velocity
H = eye(6);

% Measurement noise covariance
R = diag([gpsNoiseStd^2;
          gpsNoiseStd^2;
          gpsNoiseStd^2;
          gpsNoiseVelStd^2;
          gpsNoiseVelStd^2;
          gpsNoiseVelStd^2]);

% Process noise covariance
Q = B * (accelNoiseStd^2) * B';

% Initial state estimate
x_est = zeros(6,length(t));

x_est(:,1) = [x_true(1);
              y_true(1);
              z_true(1);
              vx_true(1);
              vy_true(1);
              vz_true(1)];

% Initial uncertainty
P = eye(6);

% Kalman filter
for k = 2:length(t)

    % IMU acceleration measurement
    a_IMU = [ax_IMU(k);
             ay_IMU(k);
             az_IMU(k)];

    % Prediction using IMU
    x_pred = F*x_est(:,k-1) + B*a_IMU;

    % Predicted uncertainty
    P_pred = F*P*F' + Q;

    % GPS position and velocity measurement
    z_GPS = [x_gps(k);
             y_gps(k);
             z_gps(k);
             vx_gps(k);
             vy_gps(k);
             vz_gps(k)];

    % Kalman gain
    K = P_pred*H'/(H*P_pred*H' + R);

    % GPS correction
    x_est(:,k) = x_pred + K*(z_GPS - H*x_pred);

    % Update uncertainty
    P = (eye(6) - K*H)*P_pred;

end

% Extract position
x_est_GPS_IMU = x_est(1,:);
y_est_GPS_IMU = x_est(2,:);
z_est_GPS_IMU = x_est(3,:);

% Extract velocity
vx_est_GPS_IMU = x_est(4,:);
vy_est_GPS_IMU = x_est(5,:);
vz_est_GPS_IMU = x_est(6,:);

% Position error
PositionError_GPS_IMU = sqrt((x_est_GPS_IMU-x_true).^2 + ...
                             (y_est_GPS_IMU-y_true).^2 + ...
                             (z_est_GPS_IMU-z_true).^2);

% Velocity error
VelError_GPS_IMU = sqrt((vx_est_GPS_IMU-vx_true).^2 + ...
                        (vy_est_GPS_IMU-vy_true).^2 + ...
                        (vz_est_GPS_IMU-vz_true).^2);

%% IMU + CAMERA

%% GPS + IMU + CAMERA

%% RTK + IMU

%% IMU + RTK + CAMERA

%% Error comparison

figure;
plot(t,PositionError_GPS, 'LineWidth',1)
hold on
plot(t,PositionError_IMU, 'LineWidth',1)
plot(t,PositionError_GPS_IMU,'Linewidth',1)
xlabel('Time [s]')
ylabel('Position Error [m]')
legend('GPS', 'IMU','GPS + IMU')
title('Position Error for Localization Architectures')

figure;
plot(t,VelError_GPS, 'LineWidth',1)
hold on
plot(t,VelError_IMU, 'LineWidth',1)
plot(t,VelError_GPS_IMU,'Linewidth',1)
xlabel('Time [s]')
ylabel('Position Error [m]')
legend('GPS', 'IMU', 'GPS + IMU')
title('Velocity Error for Localization Architectures')