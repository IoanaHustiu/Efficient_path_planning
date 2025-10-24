function xy = place2xy(p, T)
% PLACE2XY  Convierte un índice de lugar (en el grafo reducido) a [x y] en el mapa.
% p : escalar, 1..numel(T.rem_cells)
% xy: [x y], con x = columna, y = fila (sistema cartesiano Y hacia arriba)

    % Si ya has precomputado T.centr{p} = [x y], úsalo:
    if isfield(T,'centr') && numel(T.centr) >= p && ~isempty(T.centr{p})
        xy = T.centr{p};
        return;
    end

    % Si no, convierte el índice lineal del mapa completo a (fila,columna)
    lin = T.rem_cells(p);
    [H, W] = size(T.map2D); %#ok<ASGLU>
    [r, c] = ind2sub([H, W], lin);
    xy = [c, r];
end
