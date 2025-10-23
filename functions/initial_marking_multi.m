function [m0, mF, idxStart, idxGoal] = initial_marking_multi(selectedPts, B, invMap, numPlaces)
% INITIAL_MARKING_MULTI
% Genera el marcado inicial y final para VARIOS robots.
%
% INPUTS:
%   selectedPts - celda {N_r x 1} con cada elemento = [xS yS; xG yG] (0-based, cartesiano)
%                 también acepta:
%                   - matriz 2x2 (un único robot)
%                   - arreglo 2x2xN_r
%   B           - matriz binaria (1 = obstáculo, 0 = libre)
%   invMap      - vector tal que invMap(iOriginal) = iReducido (0 si obstáculo)
%   numPlaces   - número de lugares (== size(adj,1))
%
% OUTPUTS:
%   m0          - marcado inicial (numPlaces x 1) con N_r tokens repartidos en los lugares de inicio
%   mF          - marcado final   (numPlaces x 1) con N_r tokens repartidos en los lugares objetivo
%   idxStart    - (N_r x 1) índices de lugar inicial (en grafo reducido)
%   idxGoal     - (N_r x 1) índices de lugar objetivo (en grafo reducido)
%
% Ejemplo de uso:
%   selIdx      = randperm(nTotal, N_r);
%   selectedPts = robotPts(selIdx);        % robotPts como {N x 1}, cada uno [xS yS; xG yG]
%   [m0,mF,idxS,idxG] = initial_marking_multi(selectedPts, B, invMap, size(adj,1));

    % --- Normalizar entrada a celda {N_r x 1} de matrices 2x2 ---
    if iscell(selectedPts)
        robots = selectedPts;
    elseif isnumeric(selectedPts)
        if isequal(size(selectedPts), [2 2])         % un solo robot (2x2)
            robots = {double(selectedPts)};
        elseif ndims(selectedPts) == 3 && size(selectedPts,1) == 2 && size(selectedPts,2) == 2
            N_r = size(selectedPts,3);
            robots = cell(N_r,1);
            for k = 1:N_r
                robots{k} = double(selectedPts(:,:,k));
            end
        else
            error('selectedPts numérico debe ser 2x2 (1 robot) o 2x2xN_r (varios robots).');
        end
    else
        error('selectedPts debe ser celda, 2x2 o 2x2xN_r.');
    end

    N_r = numel(robots);
    if N_r == 0
        error('No hay robots en selectedPts.');
    end

    % --- Preparar salidas ---
    m0 = zeros(numPlaces, 1);
    mF = zeros(numPlaces, 1);
    idxStart = zeros(N_r, 1);
    idxGoal  = zeros(N_r, 1);

    % --- Tamaño del mapa ---
    [H, W] = size(B);

    % --- Procesar cada robot ---
    for r = 1:N_r
        R = robots{r};
        if ~isnumeric(R) || ~isequal(size(R), [2 2])
            error('Elemento %d de selectedPts no es una matriz 2x2 [xS yS; xG yG].', r);
        end
        R = double(R);

        % Coordenadas 0-based
        start_xy = R(1,:);   % [x y]
        goal_xy  = R(2,:);   % [x y]

        % --- START: a índices 1-based de MATLAB ---
        colS = start_xy(1) + 1;
        rowS = start_xy(2) + 1;
        if colS < 1 || colS > W || rowS < 1 || rowS > H
            error('Robot %d: Start fuera de rango: x=%d, y=%d', r, start_xy(1), start_xy(2));
        end
        if B(rowS, colS) == 1
            error('Robot %d: la posición inicial (%d,%d) está en un obstáculo.', r, start_xy(1), start_xy(2));
        end
        iOrgS    = sub2ind([H, W], rowS, colS);
        idxS     = invMap(iOrgS);
        if idxS == 0
            error('Robot %d: lugar inicial no encontrado en grafo reducido (invMap=0).', r);
        end

        % --- GOAL: a índices 1-based de MATLAB ---
        colG = goal_xy(1) + 1;
        rowG = goal_xy(2) + 1;
        if colG < 1 || colG > W || rowG < 1 || rowG > H
            error('Robot %d: Goal fuera de rango: x=%d, y=%d', r, goal_xy(1), goal_xy(2));
        end
        if B(rowG, colG) == 1
            error('Robot %d: la posición objetivo (%d,%d) está en un obstáculo.', r, goal_xy(1), goal_xy(2));
        end
        iOrgG    = sub2ind([H, W], rowG, colG);
        idxG     = invMap(iOrgG);
        if idxG == 0
            error('Robot %d: lugar objetivo no encontrado en grafo reducido (invMap=0).', r);
        end

        % Guardar índices
        idxStart(r) = idxS;
        idxGoal(r)  = idxG;

        % Acumular tokens en los marcados globales
        m0(idxS) = m0(idxS) + 1;
        mF(idxG) = mF(idxG) + 1;
    end

    % --- Advertencias útiles (opcionales) ---
    % Duplicidad de inicios o finales (varios robots en el mismo lugar)
    if numel(unique(idxStart)) < N_r
        warning('Hay robots con el MISMO lugar de inicio (colisiones iniciales potenciales).');
    end
    if numel(unique(idxGoal)) < N_r
        warning('Hay robots con el MISMO lugar objetivo.');
    end
end
