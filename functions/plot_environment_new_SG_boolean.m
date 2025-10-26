function env_h = plot_environment_new_SG_boolean(selectedStart, selectedFin, map2D)
% PLOT_ENVIRONMENT_NEW_SG
% Accepts starts/goals as either numeric N×2 [x y] or cell {N×1} with [x y].
% Draws free=white, obstacles=black; starts=red triangles, goals=gold diamonds.

    if nargin < 3, error('Usage: plot_environment_new_SG(S, G, map2D)'); end

    % --- normalize inputs to numeric N×2 [x y] (0-based) ---
    S = toNx2(selectedStart, 'selectedStart');
    G = toNx2(selectedFin,   'selectedFin');

    % --- ensure double map and compute normalization safely ---
    map2D = double(map2D);
    [H, W] = size(map2D);

    % match counts
    N = min(size(S,1), size(G,1));
    if N == 0, warning('No start/goal points to plot.'); env_h=[]; return; end

    % ---- figure & background (white free / black obstacle) ----
    env_h = figure('Name','Environment','NumberTitle','off'); hold on;

    cmin = min(map2D(:)); 
    cmax = max(map2D(:));
    if cmax > cmin
        map_norm = (map2D - cmin) ./ (cmax - cmin);   % in [0,1]
    else
        map_norm = zeros(size(map2D));
    end
    map_inv  = 1 - map_norm;  % obstacles black if map2D high

    surface('XData',[0 W;0 W], 'YData',[0 0;H H], 'ZData',zeros(2), ...
            'CData',map_inv, 'FaceColor','texturemap','EdgeColor','none');
    colormap(gray);
    set(gca,'YDir','normal'); axis([0 W 0 H]); axis equal tight;
    set(gcf,'Units','Normalized','OuterPosition',[0 0.5 .5 0.5]);

    % clamp to bounds (0-based coords)
    S(:,1) = min(max(S(:,1),0), W);  S(:,2) = min(max(S(:,2),0), H);
    G(:,1) = min(max(G(:,1),0), W);  G(:,2) = min(max(G(:,2),0), H);

    % ---- draw starts/goals ----
    s = 0.8; gold = [1 0.84 0];

    for k=1:size(S,1)
        c = S(k,:);
        tri = [ c(1),   c(2)+s;  c(1)-s, c(2)-s;  c(1)+s, c(2)-s ];
        patch('XData',tri(:,1),'YData',tri(:,2),'FaceColor','r','EdgeColor','k','FaceAlpha',1);
    end
    for k=1:size(G,1)
        g = G(k,:);
        dia = [ g(1),   g(2)+s;  g(1)-s, g(2);  g(1), g(2)-s;  g(1)+s, g(2) ];
        patch('XData',dia(:,1),'YData',dia(:,2),'FaceColor',gold,'EdgeColor','k','FaceAlpha',1);
    end

    xlabel('x'); ylabel('y');
    title(sprintf('Starts %i (red) & Goals %i (gold) for %d robots',N, size(G,1), N));
    grid off; hold off;
end

% ----- helper: normalize to numeric N×2 -----
function P = toNx2(inp, name)
    if iscell(inp)
        n = numel(inp);
        P = zeros(n,2);
        for i=1:n
            v = inp{i};
            if ~isnumeric(v) || numel(v)~=2
                error('%s{%d} must be numeric [x y].', name, i);
            end
            P(i,:) = v(:).';  % row [x y]
        end
    elseif isnumeric(inp)
        if size(inp,2)==2, P = double(inp);
        elseif size(inp,1)==2, P = double(inp.');  % 2×N -> N×2
        else, error('%s must be N×2 or 2×N numeric, or {N×1} cells of [x y].', name);
        end
    else
        error('%s must be numeric N×2 or cell {N×1} of [x y].', name);
    end
end
