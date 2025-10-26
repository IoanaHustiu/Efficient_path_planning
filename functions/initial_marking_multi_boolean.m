function [m0, idxStart] = initial_marking_multi_boolean(selectedStart, B, invMap, numPlaces)
% INITIAL_MARKING_MULTI_NEW_SG
% Build initial and final markings for MANY robots from separate start/goal cells.
%
% INPUTS:
%   selectedStart : cell {N×1}, each cell = [x y]  (0-based Cartesian coordinates)
%   B             : binary map (1 = obstacle, 0 = free)
%   invMap        : vector mapping full-grid linear index -> reduced-graph index (0 if obstacle)
%   numPlaces     : number of places in reduced graph
%
% OUTPUTS:
%   m0        : numPlaces×1 initial marking (token vectors)
%   idxStart      : N×1 start place indices (in reduced graph)


    % --- basic checks ----------------------------------------------------
    if nargin < 4
        error('Usage: initial_marking_multi_new(selectedStart, B, invMap, numPlaces)');
    end
    if ~iscell(selectedStart) 
        error('selectedStart must be cell arrays, each containing [x y] vectors.');
    end
    N_s = numel(selectedStart);
    N_r = min(N_s);

    [H, W] = size(B);

    % --- initialize outputs ---------------------------------------------
    m0 = zeros(numPlaces, 1);
    idxStart = zeros(N_r, 1);

    % --- process each robot ---------------------------------------------
    for r = 1:N_r
        % Extract numeric [x y]
        S = double(selectedStart{r});
        if numel(S) ~= 2 
            error('Each entry in selectedStart must be a numeric [x y].');
        end

        % --- START ---
        colS = S(1) + 1;   % 0-based → 1-based
        rowS = S(2) + 1;
        if colS < 1 || colS > W || rowS < 1 || rowS > H
            error('Robot %d: Start out of bounds (x=%d, y=%d)', r, S(1), S(2));
        end
        if B(rowS, colS) == 1
            error('Robot %d: Start (%d,%d) is on an obstacle.', r, S(1), S(2));
        end
        iOrgS = sub2ind([H, W], rowS, colS);
        idxS  = invMap(iOrgS);
        if idxS == 0
            error('Robot %d: Start not found in reduced graph (invMap=0).', r);
        end


        % --- store and update markings ---
        idxStart(r) = idxS;
        m0(idxS) = m0(idxS) + 1;
    end

    % --- warnings for duplicates ----------------------------------------
    if numel(unique(idxStart)) < N_r
        warning('Multiple robots share the same START place.');
    end
end
