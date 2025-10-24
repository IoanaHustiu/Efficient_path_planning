function [robotPts, scenTbl] = loadAllScens(baseScenName, maxNum)
% LOADALLSCENS  Reads name1.scen..nameK.scen from the current folder,
%               merges them, and removes entries with duplicate start or goal points.
%
% USAGE:
%   [robotPts, scenTbl] = loadAllScens('Paris_1_256-random.scen', 25);
%
% INPUTS:
%   baseScenName : string/char, e.g., 'Paris_1_256-random.scen'
%                  The function will look for 'Paris_1_256-random1.scen', ... '...K.scen'
%   maxNum       : positive integer, number of .scen files to read
%
% OUTPUTS:
%   robotPts     : cell array, robotPts{i} = [startX startY; goalX goalY]
%   scenTbl      : merged table with an added 'sourceFile' column
%
% REQUIREMENTS:
%   parseMovingAIScen(scenFile) must be on the MATLAB path.
%   It should return [~, T] with columns: map, startX, startY, goalX, goalY.
%
% NOTES:
% - Files are read from the CURRENT DIRECTORY (path in baseScenName is ignored).
% - Any row sharing the same start point or the same goal point with another
%   row is removed. Only the first occurrence is kept (stable order).

    % ---- Input validation
    if nargin < 2 || isempty(baseScenName) || isempty(maxNum)
        error('Usage: loadAllScens(''name.scen'', maxNum)');
    end
    if ~isscalar(maxNum) || maxNum < 1 || maxNum ~= floor(maxNum)
        error('maxNum must be a positive integer.');
    end
    baseScenName = string(baseScenName);

    % ---- Remove path and extension (we assume current folder)
    [~, baseNoPath, ext] = fileparts(baseScenName);
    if isempty(ext)
        ext = '.scen';
    end %#ok<NASGU>

    allTables = {};
    validCount = 0;
    totalRows  = 0;

    % ---- Read files name1.scen ... nameK.scen
    for k = 1:maxNum
        scenName = sprintf('%s%d.scen', baseNoPath, k);
        fid = fopen(scenName,'r');
        assert(fid>0, 'Cannot open map file: %s', scenName);
        try
            [~, t] = parseMovingAIScen(scenName);

            % Ensure 'map' is of type string for reliable merging
            if ~ismember('map', t.Properties.VariableNames)
                warning('Table from %s lacks ''map'' column. Skipped.', scenName);
                continue;
            end
            if ~isstring(t.map)
                t.map = string(t.map);
            end

            % Add the source file name as a new column
            t.sourceFile = repmat(string(scenName), height(t), 1);

            allTables{end+1} = t; %#ok<AGROW>
            validCount = validCount + 1;
            totalRows  = totalRows + height(t);
        catch ME
            warning('Error reading %s: %s (skipped).', scenName, ME.message);
        end
    end

    if isempty(allTables)
        error('No valid .scen files were found in the current folder.');
    end

    % ---- Concatenate all tables
    scenTbl = vertcat(allTables{:});

    % ---- Check required columns
    req = {'map','startX','startY','goalX','goalY'};
    missing = setdiff(req, scenTbl.Properties.VariableNames);
    if ~isempty(missing)
        error('Missing required columns: %s', strjoin(missing, ', '));
    end

    % ---- Normalize data types
    if ~isstring(scenTbl.map), scenTbl.map = string(scenTbl.map); end

    % ---------------------------------------------------------------------
    % FILTERING LOGIC:
    % Keep only the first appearance for each start and goal coordinate pair.
    % Any later row with the same (startX,startY) OR (goalX,goalY)
    % as a previous one is removed.
    % ---------------------------------------------------------------------
    n = height(scenTbl);
    keep = false(n,1);

    % Use hash maps (containers.Map) with string keys "<x>_<y>"
    seenStart = containers.Map('KeyType','char','ValueType','logical');
    seenGoal  = containers.Map('KeyType','char','ValueType','logical');

    for i = 1:n
        sx = scenTbl.startX(i); sy = scenTbl.startY(i);
        gx = scenTbl.goalX(i);  gy = scenTbl.goalY(i);

        keyS = sprintf('%d_%d', sx, sy);
        keyG = sprintf('%d_%d', gx, gy);

        if (~isKey(seenStart, keyS)) && (~isKey(seenGoal, keyG))
            keep(i) = true;
            seenStart(keyS) = true;
            seenGoal(keyG)  = true;
        else
            % Conflict: repeated start or goal → discard this row
            keep(i) = false;
        end
    end

    scenTbl = scenTbl(keep, :);

    % ---- Build robotPts output
    N = height(scenTbl);
    robotPts = cell(N,1);
    for i = 1:N
        robotPts{i} = [scenTbl.startX(i) scenTbl.startY(i); ...
                       scenTbl.goalX(i)  scenTbl.goalY(i)];
    end

    % ---- Print summary
    fprintf(['Valid files read: %d/%d. Total rows: %d. ', ...
             'Kept after removing duplicate start/goal points: %d.\n'], ...
            validCount, maxNum, totalRows, N);
end
