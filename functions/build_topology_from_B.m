function T = build_topology_from_B(B)
% BUILD_TOPOLOGY_FROM_B  Crea la topología de celdas de una rejilla a partir de B.
% 
% INPUT:
%   B : matriz binaria (H x W), 1 = obstáculo, 0 = libre
%
% OUTPUT (struct T):
%   T.adj        : matriz de adyacencia (sparse, N x N) con vecindad 4 entre celdas libres
%   T.obstacles  : índices lineales (en orden MATLAB) de celdas obstáculo (B==1)
%   T.rem_cell       : índices lineales de celdas libres (B==0)               [info extra útil]
%   T.size       : [H W]                                                 [info extra útil]
%
% Nota: La adyacencia solo conecta pares de celdas libres (no se crean aristas a/desde obstáculos).
%       Si quieres vecindad-8, ver comentario al final.

    % Validación básica
    if ~islogical(B)
        B = logical(B ~= 0);
    end
    [H, W] = size(B);
    N = H * W;

    % Máscaras e índices
    isObstacle = B;
    isFree = ~isObstacle;

    T = struct();
    T.size = [H, W];
    T.obstacles = find(isObstacle);
    T.rem_cells = find(isFree);

    % -----------------------------
    % Construcción de adyacencia 4-neighbors (N, S, E, O) SOLO entre libres
    % -----------------------------

    % Este bloque construye aristas (i<->j) sin bucles ni duplicados,
    % añadiendo solo (E) y (S) y luego simetrizamos.

    % --- Vecinos ESTE (r, c) <-> (r, c+1)
    if W > 1
        % Máscara de pares válidos (libre-libre) entre columnas adyacentes
        maskE = isFree(:, 1:W-1) & isFree(:, 2:W);

        % Coordenadas de los pares
        [rowsE, colsE1] = find(maskE);       % celdas izquierdas
        colsE2 = colsE1 + 1;                 % celdas derechas

        % Índices lineales en la matriz completa HxW
        iE = sub2ind([H, W], rowsE, colsE1);
        jE = sub2ind([H, W], rowsE, colsE2);
    else
        iE = []; jE = [];
    end

    % --- Vecinos SUR (r, c) <-> (r+1, c)
    if H > 1
        maskS = isFree(1:H-1, :) & isFree(2:H, :);

        [rowsS1, colsS] = find(maskS);       % celdas superiores
        rowsS2 = rowsS1 + 1;                 % celdas inferiores

        iS = sub2ind([H, W], rowsS1, colsS);
        jS = sub2ind([H, W], rowsS2, colsS);
    else
        iS = []; jS = [];
    end

    % Unimos aristas E y S (sin duplicar) y simetrizamos para grafo no dirigido
    I = [iE; iS; jE; jS];
    J = [jE; jS; iE; iS];

    % Creamos matriz dispersa (peso 1 por arista)
    T.adj = sparse(I, J, 1, N, N);

    % -----------------------------
    % (Opcional) Asegurar que los obstáculos no tengan aristas
    % (en principio ya se cumple por las máscaras 'isFree')
    % -----------------------------
    % T.adj(T.obstacles, :) = 0;
    % T.adj(:, T.obstacles) = 0;

end
