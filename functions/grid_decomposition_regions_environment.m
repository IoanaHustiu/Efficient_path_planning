function T = grid_decomposition_regions_environment(map_size, obs_size)
% GRID_DECOMPOSITION_REGIONS  Descompone un mapa .map en celdas regulares,
% devolviendo vértices, adyacencia, centroides y clasificación de obstáculos.
%
%   Entrada:
%     map_size : N_y, N_x, N_o
%     obs_size : longitud del lado de cada celda
%
%   Salida:
%     T.map2D      : matriz lógica [true=libre, false=obstáculo]
%     T.mapVec     : vector lógico linealizado de tamaño N
%     T.Vert       : cell(N,1) con vértices 4×2 de cada celda
%     T.adj        : sparse(N,N) adyacencia (incluye reflexivos)
%     T.centr      : cell(N,1) con centroides [x y]
%     T.obstacles  : índices lineales de celdas-obs
%     T.rem_cell   : índices lineales de celdas-libres

    %--- 1) Leer y preparar el mapa
    % % raw2D    = read_map_file(map_file);    % 0=obstáculo, 1=libre
    % % map2D    = ~logical(raw2D);            % true=libre, false=obs
    % % imagesc(map2D); colormap(gray); 
    % % set(gca, 'YDir','normal'); 
    % % axis equal tight;

    Ny = map_size(2);
    Nx = map_size(1);
    N_o = map_size(3);
    N = Nx * Ny;

    map = true(Nx, Ny);
    obstacles = randperm(N,N_o);
    map(obstacles) = false;

    map2D    = logical(map);            % true=libre, false=obs
    % imagesc(map2D); colormap(gray); 
    % set(gca, 'YDir','normal');
    % axis equal tight;

    % Vectorización y clasificación
    mapVec    = map2D(:);
    obstacles = find(~mapVec);            % celdas-obs
    rem_cell  = find(mapVec);             % celdas-libres

    %--- 2) Preparar malla de coordenadas de celda
    % Orígenes de cada celda en x e y
    [X0, Y0] = meshgrid( (0:Nx-1)*obs_size, (0:Ny-1)*obs_size );
    X0 = X0(:); Y0 = Y0(:);

    %--- 3) Prealocación de salida
    Vert  = cell(N,1);
    centr = cell(N,1);
    adj   = sparse(N,N);

    %--- 4) Rellenar vértices y centroides
    half = obs_size/2;
    for idx = 1:N
        x = X0(idx);
        y = Y0(idx);
        % Vértices [x y] en orden (bl, br, tr, tl)
        Vert{idx} = [ x,     x+obs_size, x+obs_size, x;
                      y,     y,           y+obs_size, y+obs_size ];
        centr{idx} = [ x + half, y + half ];
    end

    %--- 5) Construir adyacencia en O(N)
    [rows, cols] = ind2sub([Ny, Nx], 1:N);
    shifts = [ -1, 0; 1, 0; 0, -1; 0, 1 ];
    for k = 1:size(shifts,1)
        drow = shifts(k,1);
        dcol = shifts(k,2);
        rr = rows + drow;
        cc = cols + dcol;
        valid = rr>=1 & rr<=Ny & cc>=1 & cc<=Nx;
        src   = find(valid);
        dst   = sub2ind([Ny, Nx], rr(valid), cc(valid));
        adj   = adj + sparse(src, dst, 1, N, N);
    end
    % Reflexivos
    adj = adj + speye(N);

    %--- 6) Empaquetar resultados
    T.map2D     = map2D;
    T.mapVec    = mapVec;
    T.Vert      = Vert;
    T.centr     = centr;
    T.adj       = adj;
    T.obstacles = obstacles;
    T.rem_cell  = rem_cell;
end
