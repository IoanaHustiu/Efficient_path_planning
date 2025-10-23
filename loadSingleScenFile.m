function [robotPts, scenTbl] = loadSingleScenFile(scenFile)
% LOADSINGLESCENFILE  Load one MovingAI .scen file and remove duplicates.
%
%   [robotPts, scenTbl] = loadSingleScenFile(scenFile)
%
% INPUT:
%   scenFile - path to a MovingAI .scen file (e.g., 'Paris_1_256-random3.scen')
%
% OUTPUT:
%   robotPts - cell array {N x 1}, each = [x_start y_start; x_goal y_goal]
%   scenTbl  - table with columns: bucket, map, mapW, mapH, startX, startY, goalX, goalY, optimalLen
%
% DEPENDENCIES:
%   Requires parseMovingAIScen.m in the same path.
%
% Example:
%   [robotPts, scenTbl] = loadSingleScenFile('Paris_1_256-random3.scen');

    % --- Parse file ---
    [robotPts, scenTbl] = parseMovingAIScen(scenFile);

    % --- Remove duplicates by (map, startX, startY, goalX, goalY) ---
    [~, ia] = unique(scenTbl(:, {'map','startX','startY','goalX','goalY'}), 'rows', 'stable');
    scenTbl = scenTbl(ia, :);
    robotPts = robotPts(ia);

    % --- Optional info ---
    fprintf('Loaded %d unique start/goal pairs from "%s".\n', height(scenTbl), scenFile);
end
