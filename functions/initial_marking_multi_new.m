function [m0, mF, idxStart, idxGoal] = initial_marking_multi_new(selectedStart, selectedFin, B, invMap, numPlaces)
% INITIAL_MARKING_MULTI_NEW_SG
% Build initial and final markings for MANY robots from separate start/goal cells.
%
% INPUTS:
%   selectedStart : cell {N×1}, each cell = [x y]  (0-based Cartesian coordinates)
%   selectedFin   : cell {N×1}, each cell = [x y]  (0-based Cartesian coordinates)
%   B             : binary map (1 = obstacle, 0 = free)
%   invMap        : vector mapping full-grid linear index -> reduced-graph index (0 if obstacle)
%   numPlaces     : number of places in reduced graph
%
% OUTPUTS:
%   m0, mF        : numPlaces×1 initial/final markings (token vectors)
%   idxStart      : N×1 start place indices (in reduced graph)
%   idxGoal       : N×1 goal  place indices (in reduced graph)
%

    % --- basic checks ----------------------------------------------------
    if nargin < 5
        error('Usage: initial_marking_multi_new_SG(selectedStart, selectedFin, B, invMap, numPlaces)');
    end
    if ~iscell(selectedStart) || ~iscell(selectedFin)
        error('selectedStart and selectedFin must be cell arrays, each containing [x y] vectors.');
    end
    N_s = numel(selectedStart);
    N_g = numel(selectedFin);
    N_r = min(N_s, N_g);

    if N_s ~= N_g
        warning('Different number of starts (%d) and goals (%d). Using %d pairs.', N_s, N_g, N_r);
    end
    if N_r == 0
        error('No start/goal points provided.');
    end

    [H, W] = size(B);

    % --- initialize outputs ---------------------------------------------
    m0 = zeros(numPlaces, 1);
    mF = zeros(numPlaces, 1);
    idxStart = zeros(N_r, 1);
    idxGoal  = zeros(N_r, 1);

    % --- process each robot ---------------------------------------------
    for r = 1:N_r
        % Extract numeric [x y]
        S = double(selectedStart{r});
        G = double(selectedFin{r});
        if numel(S) ~= 2 || numel(G) ~= 2
            error('Each entry in selectedStart/selectedFin must be a numeric [x y].');
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

        % --- GOAL ---
        colG = G(1) + 1;
        rowG = G(2) + 1;
        if colG < 1 || colG > W || rowG < 1 || rowG > H
            error('Robot %d: Goal out of bounds (x=%d, y=%d)', r, G(1), G(2));
        end
        if B(rowG, colG) == 1
            error('Robot %d: Goal (%d,%d) is on an obstacle.', r, G(1), G(2));
        end
        iOrgG = sub2ind([H, W], rowG, colG);
        idxG  = invMap(iOrgG);
        if idxG == 0
            error('Robot %d: Goal not found in reduced graph (invMap=0).', r);
        end

        % --- store and update markings ---
        idxStart(r) = idxS;
        idxGoal(r)  = idxG;
        m0(idxS) = m0(idxS) + 1;
        mF(idxG) = mF(idxG) + 1;
    end

    % --- warnings for duplicates ----------------------------------------
    if numel(unique(idxStart)) < N_r
        warning('Multiple robots share the same START place.');
    end
    if numel(unique(idxGoal)) < N_r
        warning('Multiple robots share the same GOAL place.');
    end
end
