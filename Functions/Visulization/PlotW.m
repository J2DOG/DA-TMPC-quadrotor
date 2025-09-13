function PlotW(PrePoints, W_vertices, W_hat_ini_vertices)
    %% PlotW - Visualize pre-points and two convexhull (real and fitted)
    % Inputs:
    % PrePoints - N×3 matrix of trajectory points
    % W_vertices - M×3 matrix of real workspace vertices
    % W_hat_ini_vertices - K×3 matrix of fitted workspace vertices
    
    fprintf('Drawing convexhull visualization...\n');
    
    % Compute convex hull faces for both workspace surfaces
    faces_W = convhull(W_vertices);
    faces_W_hat = convhull(W_hat_ini_vertices);
    
    % Create new figure and configure properties
    fig_handle = setupFigure();
    
    % Plot trajectory with color-coded progression
    p1 = plotTrajectory(PrePoints);
    
    % Plot workspace surfaces
    [h1, h2] = plotWorkspaceSurfaces(faces_W, W_vertices, faces_W_hat, W_hat_ini_vertices);
    
    % Finalize plot appearance with proper handles
    finalizePlot(p1, h1, h2);
    
    fprintf('Convexhull visualization completed.\n');
end

function fig_handle = setupFigure()
    % Create new figure to avoid conflicts with existing plots
    fig_handle = figure;
    
    % Set up figure dimensions and properties
    width = 7.16;
    height = 3.5;
    
    set(fig_handle, 'Units', 'inches', 'Position', [1 1 width height]);
    set(fig_handle, 'PaperUnits', 'inches', ...
        'PaperPosition', [0 0 width height], ...
        'PaperSize', [width height]);
    
    % Clear any existing content
    clf(fig_handle);
    
    hold on;
    grid on;
end

function p1 = plotTrajectory(PrePoints)
    % Plot trajectory points with color progression
    x = PrePoints(:, 1);
    y = PrePoints(:, 2);
    z = PrePoints(:, 3);
    
    % Create time/index vector for color mapping
    t = 1:length(x);
    
    % Create colored trajectory line using patch
    p1 = patch(x, y, z, t, ...
        'EdgeColor', 'interp', ...
        'LineWidth', 1, ...
        'FaceColor', 'none', ...
        'DisplayName', 'Trajectory Points'); % Add display name for legend
    
    % Apply jet colormap and add colorbar
    colormap(jet);
    c = colorbar;
    c.Label.String = 'Time Progression';
end

function [h1, h2] = plotWorkspaceSurfaces(faces_W, W_vertices, faces_W_hat, W_hat_ini_vertices)
    % Plot both real and fitted workspace surfaces
    
    % Real workspace surface (gray, semi-transparent)
    h1 = trisurf(faces_W, ...
        W_vertices(:, 1), W_vertices(:, 2), W_vertices(:, 3), ...
        'FaceColor', 'g', ...
        'FaceAlpha', 0.4, ...
        'EdgeColor', 'k', ...
        'DisplayName', 'Real Workspace'); % Add display name for legend
    
    % Fitted workspace surface (green, semi-transparent)
    h2 = trisurf(faces_W_hat, ...
        W_hat_ini_vertices(:, 1), W_hat_ini_vertices(:, 2), W_hat_ini_vertices(:, 3), ...
        'FaceColor', [0.8 0.8 0.8], ...
        'FaceAlpha', 0.4, ...
        'EdgeColor', 'k', ...
        'DisplayName', 'Fitted Workspace'); % Add display name for legend
end

function finalizePlot(p1, h1, h2)
    % Set axis labels and plot properties
    xlabel('X-axis');
    ylabel('Y-axis');
    zlabel('Z-axis');
    
    % Create legend using handles to ensure correct association
    legend([p1, h1, h2], 'Location', 'best', 'FontSize', 12);
    
    % Set 3D view and enable rotation
    view(3);
    rotate3d on;
    
    % Apply equal axis scaling
    axis equal;
    
    % Properly finish the plot
    hold off;
    drawnow; % Force rendering
end