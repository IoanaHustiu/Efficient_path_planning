function [robotPts, scenTbl] = loadParis256AllScens(baseDir)
% LOADPARIS256ALLSCENS  Lee Paris_1_256-random1..25.scen, concatena y elimina duplicados.
%
% USO:
%   [robotPts, scenTbl] = loadParis256AllScens();           % en el directorio actual
%   [robotPts, scenTbl] = loadParis256AllScens('path/a/dir');
%
% Requiere en el path la función parseMovingAIScen(scenFile) proporcionada.
%
% Criterio de duplicado: mismas (map, startX, startY, goalX, goalY).

if nargin < 1 || isempty(baseDir)
    baseDir = '.';
end

allTables = {};
fileList  = cell(25,1);
rowCounts = zeros(25,1);

kValid = 0;
for k = 1:25
    scenName = sprintf('Paris_1_256-random-%d.scen', k);
    scenPath = fullfile(baseDir, scenName);
    fileList{k} = scenPath;

    if ~isfile(scenPath)
        warning('No se encuentra el archivo: %s. Se omite.', scenPath);
        continue;
    end

    try
        [~, t] = parseMovingAIScen(scenPath);
        % opcional: registrar el origen
        t.sourceFile = repmat(string(scenName), height(t), 1);
        allTables{end+1} = t; %#ok<AGROW>
        kValid = kValid + 1;
        rowCounts(k) = height(t);
    catch ME
        warning('Error leyendo %s: %s. Se omite.', scenPath, ME.message);
    end
end

if isempty(allTables)
    error('No se pudo leer ningún .scen válido en %s', baseDir);
end

% Concatenar todo
scenTbl = vertcat(allTables{:});

% Eliminar duplicados por (map, startX, startY, goalX, goalY)
[~, ia] = unique(scenTbl(:, {'map','startX','startY'}), 'rows', 'stable');
scenTbl = scenTbl(ia, :);
[~, ia] = unique(scenTbl(:, {'map','goalX','goalY'}), 'rows', 'stable');
scenTbl = scenTbl(ia, :);

% Construir robotPts con las coordenadas únicas (0-based, top-left)
N = height(scenTbl);
robotPts = cell(N,1);
for i = 1:N
    robotPts{i} = [scenTbl.startX(i) scenTbl.startY(i); scenTbl.goalX(i) scenTbl.goalY(i)];
end

% Mensaje informativo (opcional)
fprintf('Ficheros válidos: %d/25. Filas totales: %d. Únicas: %d.\n', ...
    kValid, sum(rowCounts), N);
end


