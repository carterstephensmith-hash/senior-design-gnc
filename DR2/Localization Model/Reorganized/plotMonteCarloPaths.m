function plotMonteCarloPaths(waypoints, p_true, directPaths, kfPaths, tolerance, plotTitle)

    pathBand = polybuffer(waypoints(:, 1:2), 'lines', tolerance);

    paths = {directPaths, kfPaths};
    methodNames = {'Direct RTK + IMU', 'Kalman filter'};

    % Light colors help distinguish overlapping trials
    colors = [0.55, 0.75, 0.95;
              0.95, 0.65, 0.45];

    nTrials = size(directPaths, 3);

    figure;
    tiledlayout(1, 2, 'TileSpacing', 'compact');
    sgtitle(plotTitle);

    ax = gobjects(1, 2);

    for method = 1:2

        ax(method) = nexttile;
        hold on;

        hBand = plot(pathBand, 'FaceColor', [0.70, 0.85, 0.70], 'FaceAlpha', 0.35, 'EdgeColor', 'none');

        % Each matrix column contains one trial's trajectory
        methodPaths = paths{method};
        nSamples = size(methodPaths, 2);

        east = reshape(methodPaths(1, :, :), nSamples, nTrials);
        north = reshape(methodPaths(2, :, :), nSamples, nTrials);

        hEstimates = plot(east, north, 'Color', colors(method, :), 'LineWidth', 0.5);

        % Draw the true path last so it remains visible
        hTrue = plot(p_true(1, :), p_true(2, :), 'k--', 'LineWidth', 1.8);

        xlabel('East (m)');
        ylabel('North (m)');
        title(methodNames{method});

        legend([hTrue, hBand, hEstimates(1)], {'True path', sprintf('%.0f cm horizontal path band', 100*tolerance), ...
             'Individual Monte Carlo trajectories'}, 'Location', 'best');

        axis equal;
        grid on;
        box on;

    end

    % Give both panels the same coordinate limits for fair comparison
    combinedXLimits = [ax(1).XLim, ax(2).XLim];
    combinedYLimits = [ax(1).YLim, ax(2).YLim];

    set(ax, 'XLim', [min(combinedXLimits), max(combinedXLimits)], 'YLim', [min(combinedYLimits), max(combinedYLimits)]);

end