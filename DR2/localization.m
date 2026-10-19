%% DR2, RTK + IMU localization: bias-augmented filter + Monte Carlo

%% Housekeeping
clc, clear, close all;

%% Monte Carlo settings
nSeeds = 50;

%% Rates 
imuRate = 200; % [Hz] IMU sample rate 
rtkRate = 20; % [Hz] given spec
dt = 1/imuRate;
t = 0:dt:30;
N = length(t);
rtkPeriod = round(imuRate/rtkRate);

%% Truth trajectory (independent of any estimator; computed once)
x_true = 2*cos(0.3*t);
y_true = 2*sin(0.3*t);
z_true = 2 + 0.2*sin(0.5*t);

vx_true = gradient(x_true,dt);  
vy_true = gradient(y_true,dt);  
vz_true = gradient(z_true,dt);

ax_true = gradient(vx_true,dt); 
ay_true = gradient(vy_true,dt); 
az_true = gradient(vz_true,dt);

p_true = [x_true; y_true; z_true];
v_true = [vx_true; vy_true; vz_true];
a_true = [ax_true; ay_true; az_true];

%% Sensor specs
% RTK (H-RTK ZED-F9P Ultralight)
rtkNoiseStd = 0.01; % [m], RTK position accuracy
rtkVelNoiseStd = 0.05; % [m/s], RTK velocity accuracy

accelNoiseDensity_ug_rtHz = 70; % ug/rtHz
accelNoiseStd = accelNoiseDensity_ug_rtHz*1e-6*9.81*sqrt(imuRate); % [m/s^2]

accelBiasStd = 0.002; % accelerometer bias instability [m/s^2]
biasTau = 40; % bias correlation time, [s]

%% RTK outage window (lost fix / float solution)
tDropStart = 15;  
tDropEnd = 20;
inDropout = (t >= tDropStart) & (t <= tDropEnd);

%% RTK measurement latency
rtkDelay = 0.10; % [s] assumed end-to-end RTK latency
rtkDelaySteps = round(rtkDelay/dt);

%% Requirement thresholds
posReq = 0.10; % [m]
velReq = 0.10; % [m/s]

%% Filter matrices (9-state: [p; v; b_accel])

alpha = exp(-dt/biasTau);
sigEta = accelBiasStd*sqrt(1 - alpha^2); % keeps steady-state bias std correct at this dt

F9 = [eye(3),  dt*eye(3), -0.5*dt^2*eye(3);
      zeros(3), eye(3),   -dt*eye(3);
      zeros(3), zeros(3),  alpha*eye(3)];

B9 = [0.5*dt^2*eye(3); dt*eye(3); zeros(3)]; % input is a_IMU(k)

% Process noise: white accel noise feeds p,v through the same 6x6 mapping
% as before; bias random-walk noise sits in its own block. Cross-terms
% between the accel-noise and bias blocks are neglected (approximation:
% reasonable since they're independent noise sources in the generative
% model; flagged here rather than silently assumed).

Bpv = [0.5*dt^2*eye(3); dt*eye(3)];
Qpv = Bpv*(accelNoiseStd^2)*Bpv';
Q9 = blkdiag(Qpv, sigEta^2*eye(3));

H9 = [eye(6), zeros(6,3)]; % RTK measures p,v only — never observes bias directly
R = diag([rtkNoiseStd^2*ones(3,1); rtkVelNoiseStd^2*ones(3,1)]);

StdP0 = 0.05;  
StdV0 = 0.1;
P0_9 = blkdiag(diag([StdP0^2*ones(1,3), StdV0^2*ones(1,3)]), accelBiasStd^2*eye(3));

%% Monte Carlo accumulators
sumSqPos_filt = zeros(1,N);  
sumSqVel_filt = zeros(1,N);
sumSqPos_raw  = zeros(1,N);  
sumSqVel_raw  = zeros(1,N);

rmsPosFilt_nom = zeros(nSeeds,1); 
rmsVelFilt_nom = zeros(nSeeds,1);
rmsPosFilt_drop = zeros(nSeeds,1); 
rmsVelFilt_drop = zeros(nSeeds,1);

rmsPosRaw_nom = zeros(nSeeds,1); 
rmsVelRaw_nom = zeros(nSeeds,1);
rmsPosRaw_drop = zeros(nSeeds,1); 
rmsVelRaw_drop = zeros(nSeeds,1);

pctPosFilt_nom = zeros(nSeeds,1); 
pctVelFilt_nom = zeros(nSeeds,1);
pctPosFilt_drop = zeros(nSeeds,1); 
pctVelFilt_drop = zeros(nSeeds,1);
pctPosRaw_nom = zeros(nSeeds,1); 
pctVelRaw_nom = zeros(nSeeds,1);
pctPosRaw_drop = zeros(nSeeds,1); 
pctVelRaw_drop = zeros(nSeeds,1);

consistPos_nom = zeros(nSeeds,1); % fraction of nominal samples inside the 3-sigma pos bound
consistVel_nom = zeros(nSeeds,1);

nominal = ~inDropout & (t >= 2); % exclude 2 s startup settling

rmsFn = @(e,mask) sqrt(mean(e(mask).^2));
pctFn = @(e,mask,thresh) 100*mean(e(mask) <= thresh);

for s = 1:nSeeds
    rng(s);

    % Generate accel bias
    b_accel = zeros(3,N);
    b_accel(:,1) = accelBiasStd*randn(3,1);
    for k = 2:N
        b_accel(:,k) = alpha*b_accel(:,k-1) + sigEta*randn(3,1);
    end

    % Generate IMU and RTK measurements
    a_IMU = a_true + b_accel + accelNoiseStd*randn(3,N);
    p_rtk_meas = p_true + rtkNoiseStd*randn(3,N);
    v_rtk_meas = v_true + rtkVelNoiseStd*randn(3,N);

    % Initial state
    x0_9 = [p_true(:,1); v_true(:,1); zeros(3,1)] + ...
           sqrt(diag(P0_9)).*randn(9,1);

    x_est = zeros(9,N);  x_est(:,1) = x0_9;
    P = P0_9;
    P_hist = zeros(9,N); P_hist(:,1) = diag(P0_9);

    % Bias-augmented, dual-rate, latency-compensated KF
    for k = 2:N
        x_pred = F9*x_est(:,k-1) + B9*a_IMU(:,k);
        P_pred = F9*P*F9' + Q9;

        isRtkEpoch = (mod(k-1, rtkPeriod) == 0);

        if isRtkEpoch && ~inDropout(k)
            kMeas = max(k - rtkDelaySteps, 1);
            % Latency compensation: push the stale position forward using
            % the RTK velocity measurement over the known delay, instead
            % of treating the stale sample as current.
            p_comp = p_rtk_meas(:,kMeas) + v_rtk_meas(:,kMeas)*rtkDelay;
            z = [p_comp; v_rtk_meas(:,kMeas)];

            K = P_pred*H9'/(H9*P_pred*H9' + R);
            x_est(:,k) = x_pred + K*(z - H9*x_pred);
            P = (eye(9) - K*H9)*P_pred*(eye(9) - K*H9)' + K*R*K'; % Joseph form
        else
            x_est(:,k) = x_pred;
            P = P_pred;
        end

        P_hist(:,k) = diag(P);
    end

    p_est_filt = x_est(1:3,:);
    v_est_filt = x_est(4:6,:);

    % Raw RTK-only baseline (latency-compensated the same way) 
    p_est_raw = zeros(3,N);  v_est_raw = zeros(3,N);
    p_est_raw(:,1) = p_true(:,1);  v_est_raw(:,1) = v_true(:,1);
    for k = 2:N
        isRtkEpoch = (mod(k-1, rtkPeriod) == 0);
        if isRtkEpoch && ~inDropout(k)
            kMeas = max(k - rtkDelaySteps, 1);
            p_est_raw(:,k) = p_rtk_meas(:,kMeas) + v_rtk_meas(:,kMeas)*rtkDelay;
            v_est_raw(:,k) = v_rtk_meas(:,kMeas);
        else
            p_est_raw(:,k) = p_est_raw(:,k-1);
            v_est_raw(:,k) = v_est_raw(:,k-1);
        end
    end

    % Errors
    PosErr_filt = vecnorm(p_est_filt - p_true, 2, 1);
    VelErr_filt = vecnorm(v_est_filt - v_true, 2, 1);
    PosErr_raw = vecnorm(p_est_raw - p_true, 2, 1);
    VelErr_raw = vecnorm(v_est_raw - v_true, 2, 1);

    sumSqPos_filt = sumSqPos_filt + PosErr_filt.^2;
    sumSqVel_filt = sumSqVel_filt + VelErr_filt.^2;
    sumSqPos_raw = sumSqPos_raw + PosErr_raw.^2;
    sumSqVel_raw = sumSqVel_raw + VelErr_raw.^2;

    rmsPosFilt_nom(s) = rmsFn(PosErr_filt, nominal);
    rmsVelFilt_nom(s) = rmsFn(VelErr_filt, nominal);
    rmsPosFilt_drop(s) = rmsFn(PosErr_filt, inDropout);
    rmsVelFilt_drop(s) = rmsFn(VelErr_filt, inDropout);
    rmsPosRaw_nom(s) = rmsFn(PosErr_raw, nominal);
    rmsVelRaw_nom(s) = rmsFn(VelErr_raw, nominal);
    rmsPosRaw_drop(s) = rmsFn(PosErr_raw, inDropout);
    rmsVelRaw_drop(s) = rmsFn(VelErr_raw, inDropout);

    pctPosFilt_nom(s) = pctFn(PosErr_filt, nominal, posReq);
    pctVelFilt_nom(s) = pctFn(VelErr_filt, nominal, velReq);
    pctPosFilt_drop(s) = pctFn(PosErr_filt, inDropout, posReq);
    pctVelFilt_drop(s) = pctFn(VelErr_filt, inDropout, velReq);
    pctPosRaw_nom(s) = pctFn(PosErr_raw, nominal, posReq);
    pctVelRaw_nom(s) = pctFn(VelErr_raw, nominal, velReq);
    pctPosRaw_drop(s) = pctFn(PosErr_raw, inDropout, posReq);
    pctVelRaw_drop(s) = pctFn(VelErr_raw, inDropout, velReq);

    posSigma3 = 3*sqrt(sum(P_hist(1:3,:),1));
    velSigma3 = 3*sqrt(sum(P_hist(4:6,:),1));
    consistPos_nom(s) = 100*mean(PosErr_filt(nominal) <= posSigma3(nominal));
    consistVel_nom(s) = 100*mean(VelErr_filt(nominal) <= velSigma3(nominal));

    if s == 1
        % Keep seed-1 time histories for the plots below (matches v1's
        % single-run plots, but from the improved filter).
        PosErr_filt_s1 = PosErr_filt; 
        VelErr_filt_s1 = VelErr_filt;
        PosErr_raw_s1 = PosErr_raw;  
        VelErr_raw_s1 = VelErr_raw;
        posSigma3_s1 = posSigma3;
        velSigma3_s1 = velSigma3;
    end
end

RMSE_Pos_filt_t = sqrt(sumSqPos_filt/nSeeds); % expected error vs time, across seeds
RMSE_Vel_filt_t = sqrt(sumSqVel_filt/nSeeds);
RMSE_Pos_raw_t = sqrt(sumSqPos_raw/nSeeds);
RMSE_Vel_raw_t = sqrt(sumSqVel_raw/nSeeds);

%% Plots

% Expected error vs time (RMSE across all seeds), dropout shaded
figure('Name','Expected error vs time (Monte Carlo RMSE)');
tl = tiledlayout(2,1,'TileSpacing','compact');

ax1 = nexttile(tl); hold(ax1,'on')
patch(ax1,[tDropStart tDropEnd tDropEnd tDropStart],[1e-4 1e-4 10 10], ...
      [1 0.9 0.9],'EdgeColor','none','FaceAlpha',0.6)
plot(ax1,t,RMSE_Pos_raw_t,'Color',[0.6 0.6 0.6],'LineWidth',1)
plot(ax1,t,RMSE_Pos_filt_t,'b','LineWidth',1.3)
yline(ax1,posReq,'r--','LineWidth',1.2)
set(ax1,'YScale','log'); grid(ax1,'on'); ylim(ax1,[1e-4 10])
xlabel(ax1,'Time [s]'); ylabel(ax1,'Position RMSE [m]')
legend(ax1,{'RTK outage','Raw RTK','Filtered (9-state RTK+IMU)','10 cm requirement'}, ...
       'Location','northoutside','Orientation','horizontal')

ax2 = nexttile(tl); hold(ax2,'on')
patch(ax2,[tDropStart tDropEnd tDropEnd tDropStart],[1e-4 1e-4 10 10], ...
      [1 0.9 0.9],'EdgeColor','none','FaceAlpha',0.6)
plot(ax2,t,RMSE_Vel_raw_t,'Color',[0.6 0.6 0.6],'LineWidth',1)
plot(ax2,t,RMSE_Vel_filt_t,'b','LineWidth',1.3)
yline(ax2,velReq,'r--','LineWidth',1.2)
set(ax2,'YScale','log'); grid(ax2,'on'); ylim(ax2,[1e-4 10])
xlabel(ax2,'Time [s]'); ylabel(ax2,'Velocity RMSE [m/s]')

%%  Summary tables

RmsTable = table( ...
    [mean(rmsPosFilt_nom); mean(rmsPosFilt_drop); mean(rmsPosRaw_nom); mean(rmsPosRaw_drop)], ...
    [std(rmsPosFilt_nom);  std(rmsPosFilt_drop);  std(rmsPosRaw_nom);  std(rmsPosRaw_drop)], ...
    [mean(rmsVelFilt_nom); mean(rmsVelFilt_drop); mean(rmsVelRaw_nom); mean(rmsVelRaw_drop)], ...
    [std(rmsVelFilt_nom);  std(rmsVelFilt_drop);  std(rmsVelRaw_nom);  std(rmsVelRaw_drop)], ...
    [mean(pctPosFilt_nom); mean(pctPosFilt_drop); mean(pctPosRaw_nom); mean(pctPosRaw_drop)], ...
    [mean(pctVelFilt_nom); mean(pctVelFilt_drop); mean(pctVelRaw_nom); mean(pctVelRaw_drop)], ...
    'RowNames', {'Filtered - nominal','Filtered - dropout','Raw RTK - nominal','Raw RTK - dropout'}, ...
    'VariableNames', {'PosRMS_mean_m','PosRMS_std_m','VelRMS_mean_mps','VelRMS_std_mps', ...
                      'PctMeetingPosReq_mean','PctMeetingVelReq_mean'});

fprintf('\n=== DR2 results, %d-seed Monte Carlo ===\n', nSeeds);
disp(RmsTable)

ConsistTable = table(mean(consistPos_nom), std(consistPos_nom), ...
                      mean(consistVel_nom), std(consistVel_nom), ...
    'VariableNames', {'PctInside3SigmaPos_mean','PctInside3SigmaPos_std', ...
                      'PctInside3SigmaVel_mean','PctInside3SigmaVel_std'});
%% Path (DR2 inset lawnmower and kinematics, matches DR2/coveragePlanning.m)

cov.areaWidth  = 5; % [m]
cov.areaLength = 5; % [m]
cov.speed = 0.40; % [m/s] nominal cruise
cov.accel = 0.50; % [m/s^2] commanded accel/decel on each leg (TBD)
cov.cornerHold = 1.0; % [s] stop-and-settle at each waypoint (TBD)
cov.hFOV = 60; % [deg]
cov.imageResolution = [1920, 1080]; % [px] (TBD: selected camera)
cov.altitude = 2; % [m]
cov.overlap = 0.20;
cov.markerSize = 0.08; % [m]

% 3-sigma margins: 5 cm position error + 1 deg attitude jitter, plus tilt while accelerating
cov.crossTrackMargin = 3*0.05 + cov.altitude*tand(3*1.0);
cov.alongTrackMargin = cov.crossTrackMargin + cov.altitude*cov.accel/9.81;

waypoints = generateInsetLawnmower(cov);
traj = buildTrajectory(waypoints, cov.speed, cov.accel, cov.cornerHold, 0.01);
searchTime_path = traj.t(end);

% Resample the commanded trajectory onto the IMU time grid; truth velocity and
% acceleration come from the trapezoidal profile, not numerical differentiation
t_path = 0:dt:searchTime_path;
N_path = length(t_path);
p_true_path = [interp1(traj.t, traj.pos, t_path)'; cov.altitude*ones(1,N_path)];
v_true_path = [interp1(traj.t, traj.vel, t_path)'; zeros(1,N_path)];
a_true_path = [interp1(traj.t, traj.acc, t_path)'; zeros(1,N_path)];

distAlongPath = [0, cumsum(vecnorm(diff(p_true_path,1,2),2,1))];

tDropStart_path = 0.40*searchTime_path;
tDropEnd_path = tDropStart_path + 5;
inDropout_path = (t_path >= tDropStart_path) & (t_path <= tDropEnd_path);
nominal_path = ~inDropout_path & (t_path >= 2);

nSeedsPath = 50;

sumSqPos_filt_path = zeros(1,N_path);  sumSqVel_filt_path = zeros(1,N_path);
sumSqPos_raw_path  = zeros(1,N_path);  sumSqVel_raw_path  = zeros(1,N_path);

rmsPosFilt_nom_p  = zeros(nSeedsPath,1); rmsVelFilt_nom_p  = zeros(nSeedsPath,1);
rmsPosFilt_drop_p = zeros(nSeedsPath,1); rmsVelFilt_drop_p = zeros(nSeedsPath,1);
rmsPosRaw_nom_p   = zeros(nSeedsPath,1); rmsVelRaw_nom_p   = zeros(nSeedsPath,1);
rmsPosRaw_drop_p  = zeros(nSeedsPath,1); rmsVelRaw_drop_p  = zeros(nSeedsPath,1);

pctPosFilt_nom_p  = zeros(nSeedsPath,1); pctVelFilt_nom_p  = zeros(nSeedsPath,1);
pctPosFilt_drop_p = zeros(nSeedsPath,1); pctVelFilt_drop_p = zeros(nSeedsPath,1);
pctPosRaw_nom_p   = zeros(nSeedsPath,1); pctVelRaw_nom_p   = zeros(nSeedsPath,1);
pctPosRaw_drop_p  = zeros(nSeedsPath,1); pctVelRaw_drop_p  = zeros(nSeedsPath,1);

for s = 1:nSeedsPath
    rng(1000+s);

    b_accel = zeros(3,N_path);
    b_accel(:,1) = accelBiasStd*randn(3,1);
    for k = 2:N_path
        b_accel(:,k) = alpha*b_accel(:,k-1) + sigEta*randn(3,1);
    end

    a_IMU_path = a_true_path + b_accel + accelNoiseStd*randn(3,N_path);
    p_rtk_path = p_true_path + rtkNoiseStd*randn(3,N_path);
    v_rtk_path = v_true_path + rtkVelNoiseStd*randn(3,N_path);

    x0 = [p_true_path(:,1); v_true_path(:,1); zeros(3,1)] + sqrt(diag(P0_9)).*randn(9,1);
    x_est = zeros(9,N_path);  x_est(:,1) = x0;
    P = P0_9;

    p_est_raw_path = zeros(3,N_path);  v_est_raw_path = zeros(3,N_path);
    p_est_raw_path(:,1) = p_true_path(:,1);  v_est_raw_path(:,1) = v_true_path(:,1);

    for k = 2:N_path
        x_pred = F9*x_est(:,k-1) + B9*a_IMU_path(:,k);
        P_pred = F9*P*F9' + Q9;

        isRtkEpoch = (mod(k-1, rtkPeriod) == 0);
        if isRtkEpoch && ~inDropout_path(k)
            kMeas = max(k - rtkDelaySteps, 1);
            p_comp = p_rtk_path(:,kMeas) + v_rtk_path(:,kMeas)*rtkDelay;
            z = [p_comp; v_rtk_path(:,kMeas)];

            K = P_pred*H9'/(H9*P_pred*H9' + R);
            x_est(:,k) = x_pred + K*(z - H9*x_pred);
            P = (eye(9) - K*H9)*P_pred*(eye(9) - K*H9)' + K*R*K';

            p_est_raw_path(:,k) = p_comp;
            v_est_raw_path(:,k) = v_rtk_path(:,kMeas);
        else
            x_est(:,k) = x_pred;
            P = P_pred;

            p_est_raw_path(:,k) = p_est_raw_path(:,k-1);
            v_est_raw_path(:,k) = v_est_raw_path(:,k-1);
        end
    end

    PosErr_filt_p = vecnorm(x_est(1:3,:) - p_true_path, 2, 1);
    VelErr_filt_p = vecnorm(x_est(4:6,:) - v_true_path, 2, 1);
    PosErr_raw_p  = vecnorm(p_est_raw_path - p_true_path, 2, 1);
    VelErr_raw_p  = vecnorm(v_est_raw_path - v_true_path, 2, 1);

    sumSqPos_filt_path = sumSqPos_filt_path + PosErr_filt_p.^2;
    sumSqVel_filt_path = sumSqVel_filt_path + VelErr_filt_p.^2;
    sumSqPos_raw_path  = sumSqPos_raw_path  + PosErr_raw_p.^2;
    sumSqVel_raw_path  = sumSqVel_raw_path  + VelErr_raw_p.^2;

    rmsPosFilt_nom_p(s)  = rmsFn(PosErr_filt_p, nominal_path);
    rmsVelFilt_nom_p(s)  = rmsFn(VelErr_filt_p, nominal_path);
    rmsPosFilt_drop_p(s) = rmsFn(PosErr_filt_p, inDropout_path);
    rmsVelFilt_drop_p(s) = rmsFn(VelErr_filt_p, inDropout_path);
    rmsPosRaw_nom_p(s)   = rmsFn(PosErr_raw_p, nominal_path);
    rmsVelRaw_nom_p(s)   = rmsFn(VelErr_raw_p, nominal_path);
    rmsPosRaw_drop_p(s)  = rmsFn(PosErr_raw_p, inDropout_path);
    rmsVelRaw_drop_p(s)  = rmsFn(VelErr_raw_p, inDropout_path);

    pctPosFilt_nom_p(s)  = pctFn(PosErr_filt_p, nominal_path, posReq);
    pctVelFilt_nom_p(s)  = pctFn(VelErr_filt_p, nominal_path, velReq);
    pctPosFilt_drop_p(s) = pctFn(PosErr_filt_p, inDropout_path, posReq);
    pctVelFilt_drop_p(s) = pctFn(VelErr_filt_p, inDropout_path, velReq);
    pctPosRaw_nom_p(s)   = pctFn(PosErr_raw_p, nominal_path, posReq);
    pctVelRaw_nom_p(s)   = pctFn(VelErr_raw_p, nominal_path, velReq);
    pctPosRaw_drop_p(s)  = pctFn(PosErr_raw_p, inDropout_path, posReq);
    pctVelRaw_drop_p(s)  = pctFn(VelErr_raw_p, inDropout_path, velReq);

    if s == 1
        % Keep seed-1's actual estimated trajectories for the spatial plot
        p_est_filt_path_s1 = x_est(1:3,:);
        p_est_raw_path_s1  = p_est_raw_path;
    end
end

RMSE_Pos_filt_path = sqrt(sumSqPos_filt_path/nSeedsPath);
RMSE_Vel_filt_path = sqrt(sumSqVel_filt_path/nSeedsPath);
RMSE_Pos_raw_path  = sqrt(sumSqPos_raw_path/nSeedsPath);
RMSE_Vel_raw_path  = sqrt(sumSqVel_raw_path/nSeedsPath);


figure('Name','Lawnmower path — spatial trajectories, seed 1');
hold on
plot(p_true_path(1,:), p_true_path(2,:), 'g', 'LineWidth', 2)
plot(p_est_raw_path_s1(1,:), p_est_raw_path_s1(2,:), 'm', 'LineWidth', 1)
plot(p_est_filt_path_s1(1,:), p_est_filt_path_s1(2,:), 'b', 'LineWidth', 0.5)
plot(p_true_path(1,inDropout_path), p_true_path(2,inDropout_path), 'r--', 'LineWidth', 3)
axis equal; grid on; box on
xlabel('X [m]'); ylabel('Y [m]')
title(sprintf('RTK outage at t = %.0f-%.0f s (red segment)', tDropStart_path, tDropEnd_path))
legend({'True path','Raw RTK estimate','Filtered estimate','RTK outage segment'}, ...
       'Location','bestoutside')


% DR2 inset lawnmower (same layout as DR2/coveragePlanning.m): tracks inset so
% the footprint, shrunk by half a marker and the 3-sigma margins, still reaches
% every possible marker centre
function waypoints = generateInsetLawnmower(cov)
    vFOV = 2*atand((cov.imageResolution(2)/cov.imageResolution(1))*tand(cov.hFOV/2));
    half = cov.markerSize/2;
    halfWidth = cov.altitude*tand(cov.hFOV/2) - half - cov.crossTrackMargin;
    halfLength = cov.altitude*tand(vFOV/2) - half - cov.alongTrackMargin;

    span = cov.areaWidth - 2*half;
    trackSpacing = 2*halfWidth*(1 - cov.overlap);
    nTracks = max(1, ceil((span - 2*halfWidth)/trackSpacing) + 1);

    if nTracks == 1
        xlines = cov.areaWidth/2;
    else
        xlines = linspace(half + halfWidth, cov.areaWidth - half - halfWidth, nTracks);
    end

    yLow = min(half + halfLength, cov.areaLength/2);
    yHigh = max(cov.areaLength - half - halfLength, cov.areaLength/2);

    waypoints = [];
    for i = 1:nTracks
        if mod(i,2) == 1
            ys = [yLow, yHigh];
        else
            ys = [yHigh, yLow];
        end
        waypoints = [waypoints; xlines(i), ys(1); xlines(i), ys(2)];
    end
end

% Each leg starts and ends at rest with a trapezoidal speed profile, followed by
% a settle hold at the waypoint (same as DR2/coveragePlanning.m)
function traj = buildTrajectory(waypoints, speed, accel, cornerHold, dt)
    t = 0; pos = waypoints(1,:); vel = [0 0]; acc = [0 0];
    tNow = 0;

    for i = 1:size(waypoints,1)-1
        startPoint = waypoints(i,:);
        endPoint = waypoints(i+1,:);
        legLength = norm(endPoint - startPoint);
        if legLength < 1e-10, continue; end
        direction = (endPoint - startPoint)/legLength;

        [tau, s, sd, sdd] = trapezoidProfile(legLength, speed, accel, dt);
        tau = tau(2:end); s = s(2:end); sd = sd(2:end); sdd = sdd(2:end);

        t = [t; tNow + tau];
        pos = [pos; startPoint + s.*direction];
        vel = [vel; sd.*direction];
        acc = [acc; sdd.*direction];
        tNow = t(end);

        nHold = round(cornerHold/dt);
        if nHold > 0
            t = [t; tNow + (1:nHold)'*dt];
            pos = [pos; repmat(endPoint, nHold, 1)];
            vel = [vel; zeros(nHold, 2)];
            acc = [acc; zeros(nHold, 2)];
            tNow = t(end);
        end
    end

    traj.t = t; traj.pos = pos; traj.vel = vel; traj.acc = acc;
end

function [tau, s, sd, sdd] = trapezoidProfile(legLength, speed, accel, dt)
    tAccel = speed/accel;
    dAccel = 0.5*accel*tAccel^2;

    if 2*dAccel > legLength      % triangular: never reaches cruise speed
        tAccel = sqrt(legLength/accel);
        speed = accel*tAccel;
        tCruise = 0;
    else
        tCruise = (legLength - 2*dAccel)/speed;
    end

    T = 2*tAccel + tCruise;
    tau = (0:dt:T)';
    if T - tau(end) > 1e-9, tau = [tau; T]; end

    s = zeros(size(tau)); sd = s; sdd = s;
    phase1 = tau <= tAccel;
    phase2 = tau > tAccel & tau <= tAccel + tCruise;
    phase3 = tau > tAccel + tCruise;

    s(phase1) = 0.5*accel*tau(phase1).^2;
    sd(phase1) = accel*tau(phase1);
    sdd(phase1) = accel;

    s(phase2) = 0.5*accel*tAccel^2 + speed*(tau(phase2) - tAccel);
    sd(phase2) = speed;

    td = tau(phase3) - tAccel - tCruise;
    s(phase3) = 0.5*accel*tAccel^2 + speed*tCruise + speed*td - 0.5*accel*td.^2;
    sd(phase3) = speed - accel*td;
    sdd(phase3) = -accel;
end
