function B = load_map(mapFile)
% LOAD_MAP_AND_SCEN
%   Reads a MovingAI .map file into a binary matrix (1=obstacle, 0=free),
%   parses the .scen file (start/goal positions),
%   and plots the map using a Cartesian coordinate system (0,0 bottom-left).
%
%   Inputs:
%       mapFile  - path to .map file
%       scenFile - path to .scen file
%
%   Outputs:
%       B         - binary occupancy matrix (1=obstacle, 0=free)
%       robotPts  - cell array, robotPts{i} = [sx sy; gx gy] (Cartesian coords)
%       rawScen   - table with parsed .scen columns
%
% Example:
%   [B, robotPts] = load_map_and_scen('Paris_1_512.map', 'Paris_1_512.map.scen');

if nargin < 3, showK = 0; end

% ---------- Read .map ----------
B = parseMovingAIMap(mapFile);

end

% =======================================================
function B = parseMovingAIMap(mapFile)
% Reads a MovingAI .map file and returns a binary matrix B
% '@' = obstacle (1)
% others = free (0)

fid = fopen(mapFile,'r');
assert(fid>0, 'Cannot open map file: %s', mapFile);
cleanup = onCleanup(@() fclose(fid));

rows = {};
headerDone = false;
while ~feof(fid)
    line = fgetl(fid);
    if ~ischar(line), break; end
    if headerDone
        rows{end+1} = line; %#ok<AGROW>
    elseif strcmpi(strtrim(line), 'map')
        headerDone = true;
    end
end

assert(~isempty(rows), 'No map data found in file: %s', mapFile);

H = numel(rows);
W = numel(rows{1});
B = zeros(H, W, 'uint8');

for r = 1:H
    rowStr = rows{r};
    if numel(rowStr) ~= W
        error('Inconsistent row width at row %d', r);
    end
    B(r, :) = uint8(rowStr == '@' | rowStr == 'T');
end
end


