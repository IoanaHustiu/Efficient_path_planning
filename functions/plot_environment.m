function env_h = plot_environment(initial_regions, final_regions, map2D, env_limit, cells)
% PLOT_ENVIRONMENT Dibuja el entorno con celdas, obstáculos, inicios y metas
%
%   env_h           : handle de la figura
%   initial_regions : índices de celdas de inicio (lineales)
%   final_regions   : índices de celdas de meta (lineales)
%   obstacles       : índices de celdas-obstáculo (lineales)
%   env_limit       : tamaño de la rejilla (número de celdas por lado)
%   cells           : cell array con vértices 4×2 de cada celda

env_h = figure();
hold on;
imagesc([0 size(map2D,2)],[0 size(map2D,1)],map2D); colormap(gray);

% Asegura orientación tipo imagen y límites exactos
set(gca, 'YDir','normal');
axis equal tight;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0.5 .5 0.5]);


% Función auxiliar para regiones coloreadas
    function plot_regions(list, color)
        for k = 1:numel(list)
            idx = list(k);
            poly = cells{idx};
            fill(poly(1,:), poly(2,:), color, ...
                'FaceAlpha', 0.9, 'EdgeColor', 'k');
        end
    end
% Obstáculos en negro
% Inicios en rojo con etiqueta
%plot_regions(initial_regions, [1 0 0]);

% Metas en azul con etiqueta
%plot_regions(final_regions, [0 0 1]);

% 2) Iconos de inicio: flechas (triángulos rojos)
for k = 1:numel(initial_regions)
    idx = initial_regions(k);
    c = mean(cells{idx},2)';  % centro
    s = 1;
    arrow = [c(1),   c(2)+s;        % punta superior
        c(1)-s, c(2)-s;        % esquina inferior izquierda
        c(1)+s, c(2)-s];       % esquina inferior derecha
    patch(arrow(:,1), arrow(:,2), 'r', 'EdgeColor','k', 'FaceAlpha',1);
end

% 3) Iconos de meta: diamantes dorados
goldColor = [1, 0.84, 0];
for k = 1:numel(final_regions)
    idx = final_regions(k);
    c = mean(cells{idx},2)';
    s = 1;
    diamond = [c(1),   c(2)+s;       % arriba
        c(1)-s, c(2);        % izquierda
        c(1),   c(2)-s;       % abajo
        c(1)+s, c(2)];       % derecha
    patch(diamond(:,1), diamond(:,2), goldColor, 'EdgeColor','k', 'FaceAlpha',1);
    
    % centr=mean(cells{final_regions(k)},2)';
    % text(centr(1),centr(2),sprintf('g_{%d}',k),'HorizontalAlignment','center','Color','black','FontSize',7,'FontName','Times New Roman');
end
hold off;
end


