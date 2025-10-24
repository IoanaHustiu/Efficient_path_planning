function env_h = plot_environment_new(selectedPts, map2D, T)
% PLOT_ENVIRONMENT_NEW
% Draws the environment with free space = white and obstacles = black.
% Robot starts are red triangles; goals are gold diamonds.
%
%   env_h = plot_environment_new(selectedPts, map2D, T)
%
%   Inputs:
%     selectedPts : cell array, each cell is 2x2 [x_start y_start; x_goal y_goal]
%     map2D       : numeric/logical matrix (already flipped to Cartesian)
%     T           : struct (optional, kept for compatibility)
%
%   Output:
%     env_h : figure handle
%
%   Author: GPT-5 (2025) – revised for white free cells, black obstacles
% -------------------------------------------------------------------------

    % --- Input checks ----------------------------------------------------
    if nargin < 3
        error('Usage: plot_environment_new(selectedPts, map2D, T)');
    end
    if ~isnumeric(map2D)
        map2D = double(map2D);
    end
    [H, W] = size(map2D);

    % --- Figure setup ----------------------------------------------------
    env_h = figure('Name','Environment','NumberTitle','off');
    hold on;

    % --- Normalize and invert map ---------------------------------------
    % map2D = 0 (free) → white; map2D = 1 (obstacle) → black
    cmin = min(map2D(:));
    cmax = max(map2D(:));
    if cmax > cmin
        map_norm = (map2D - cmin) / (cmax - cmin);
    else
        map_norm = zeros(size(map2D));
    end
    map_inv = 1 - map_norm;   % invert colors: obstacles black, free white

    % --- Render background (no imagesc/image to avoid shadowing) ---------
    xSurf = [0, W; 0, W];
    ySurf = [0, 0; H, H];
    zSurf = zeros(2, 2);
    surface('XData', xSurf, 'YData', ySurf, 'ZData', zSurf, ...
            'CData', map_inv, ...
            'FaceColor', 'texturemap', ...
            'EdgeColor', 'none');

    colormap(gray);           % gray → 0 black, 1 white
    set(gca, 'YDir', 'normal');  % Cartesian coordinates
    axis([0 W 0 H]);
    axis equal tight;
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0.5 .5 0.5]);

    % --- Extract and validate start/goal points -------------------------
    n = numel(selectedPts);
    starts = nan(n,2);
    goals  = nan(n,2);
    for k = 1:n
        pts = selectedPts{k};
        if isnumeric(pts) && size(pts,1) >= 2 && size(pts,2) == 2
            starts(k,:) = pts(1,:);
            goals(k,:)  = pts(2,:);
        end
    end
    valid = all(isfinite(starts),2) & all(isfinite(goals),2);
    starts = starts(valid,:);
    goals  = goals(valid,:);

    starts(:,1) = min(max(starts(:,1), 0), W);
    starts(:,2) = min(max(starts(:,2), 0), H);
    goals(:,1)  = min(max(goals(:,1), 0),  W);
    goals(:,2)  = min(max(goals(:,2), 0),  H);

    % --- Plot starts/goals ----------------------------------------------
    s = 0.8;                 % icon half-size
    gold = [1, 0.84, 0];

    % Red triangles (starts)
    for k = 1:size(starts,1)
        c = starts(k,:);
        tri = [ c(1),   c(2)+s;    % up
                c(1)-s, c(2)-s;    % left
                c(1)+s, c(2)-s ];  % right
        patch('XData', tri(:,1), 'YData', tri(:,2), ...
              'FaceColor', 'r', 'EdgeColor', 'k', 'FaceAlpha', 1);
    end

    % Gold diamonds (goals)
    for k = 1:size(goals,1)
        g = goals(k,:);
        dia = [ g(1),   g(2)+s;    % up
                g(1)-s, g(2);      % left
                g(1),   g(2)-s;    % down
                g(1)+s, g(2) ];    % right
        patch('XData', dia(:,1), 'YData', dia(:,2), ...
              'FaceColor', gold, 'EdgeColor', 'k', 'FaceAlpha', 1);
    end

    % --- Labels & finishing touches -------------------------------------
    xlabel('x');
    ylabel('y');
    title(sprintf('Start (red) & Goal (gold) positions for %d robots', numel(selectedPts)));
    grid off;
    hold off;
end
