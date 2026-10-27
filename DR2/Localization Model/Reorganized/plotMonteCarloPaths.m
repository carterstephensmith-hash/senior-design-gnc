function plotMonteCarloPaths(waypoints, p_true, directPaths, kfPaths, ...
    tolerance, plotTitle)

    %% Setup

    pathBand = polybuffer(waypoints(:,1:2), 'lines', tolerance);

    paths = {directPaths, kfPaths};
    methodNames = {'Direct RTK + IMU', 'Kalman Filter'};

    % Colors
    colors = [0.25, 0.55, 0.85;    % Direct RTK + IMU
              0.90, 0.40, 0.20];   % Kalman filter

    bandColor = [0.70, 0.85, 0.70];
    trueColor = [0.10, 0.10, 0.10];

    nTrials = size(directPaths,3);

    %% Plot Monte Carlo trajectories

    figure('Name',plotTitle,'NumberTitle','off','Color','w');
    
    tl = tiledlayout(1,2, ...
        'TileSpacing','compact', ...
        'Padding','compact');

    sgtitle(tl,plotTitle, ...
        'FontSize',14, ...
        'FontWeight','bold');

    ax = gobjects(1,2);

    for method = 1:2

        ax(method) = nexttile;
        hold(ax(method),'on');

        % Allowable horizontal path corridor
        hBand = plot(ax(method),pathBand, ...
            'FaceColor',bandColor, ...
            'FaceAlpha',0.35, ...
            'EdgeColor','none');

        % Extract individual Monte Carlo trajectories
        methodPaths = paths{method};
        nSamples = size(methodPaths,2);

        east = reshape(methodPaths(1,:,:),nSamples,nTrials);
        north = reshape(methodPaths(2,:,:),nSamples,nTrials);

        % Plot all trial trajectories
        hEstimates = plot(ax(method),east,north, ...
            'Color',colors(method,:), ...
            'LineWidth',0.7);

        % True trajectory
        hTrue = plot(ax(method),p_true(1,:),p_true(2,:), ...
            '--', ...
            'Color',trueColor, ...
            'LineWidth',2);

        % Labels and formatting
        xlabel(ax(method),'East (m)');
        ylabel(ax(method),'North (m)');
        title(ax(method),methodNames{method}, ...
            'FontSize',12, ...
            'FontWeight','bold');

        grid(ax(method),'on');
        box(ax(method),'on');
        axis(ax(method),'equal');

        set(ax(method), ...
            'FontSize',11, ...
            'LineWidth',0.8, ...
            'Layer','top');

        legend(ax(method), ...
            [hTrue,hBand,hEstimates(1)], ...
            {'True path', ...
             sprintf('Path corridor (\\pm %.0f cm)',100*tolerance), ...
             'Monte Carlo trajectories'}, ...
            'Location','best', ...
            'FontSize',9);

    end

    %% Match coordinate limits across both panels

    combinedXLimits = [ax(1).XLim, ax(2).XLim];
    combinedYLimits = [ax(1).YLim, ax(2).YLim];

    commonXLimits = [min(combinedXLimits),max(combinedXLimits)];
    commonYLimits = [min(combinedYLimits),max(combinedYLimits)];

    set(ax, ...
        'XLim',commonXLimits, ...
        'YLim',commonYLimits);

end
