function [robotPts, scenTbl] = parseMovingAIScen(scenFile)
% Reads a MovingAI .scen file.
% Format: bucket mapName mapW mapH startX startY goalX goalY optimalLength
% Output: robotPts{i} = [sx sy; gx gy] (0-based top-left coords)

fid = fopen(scenFile,'r');
assert(fid>0, 'Cannot open scen file: %s', scenFile);
cleanup = onCleanup(@() fclose(fid));

firstLine = fgetl(fid);
if ~ischar(firstLine)
    error('Empty .scen file: %s', scenFile);
end
if startsWith(strtrim(lower(firstLine)), 'version')
    % skip version line
else
    fseek(fid, 0, 'bof'); % rewind if not version header
end

C = textscan(fid, '%d %s %d %d %d %d %d %d %f', ...
    'Delimiter', {' ','\t'}, 'MultipleDelimsAsOne', true);

scenTbl = table(C{1}, C{2}, C{3}, C{4}, C{5}, C{6}, C{7}, C{8}, C{9}, ...
    'VariableNames', {'bucket','map','mapW','mapH','startX','startY','goalX','goalY','optimalLen'});

N = height(scenTbl);
robotPts = cell(N,1);
for i = 1:N
    robotPts{i} = [scenTbl.startX(i) scenTbl.startY(i); scenTbl.goalX(i) scenTbl.goalY(i)];
end
end