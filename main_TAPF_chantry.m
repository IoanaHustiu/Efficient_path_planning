close all;
clear;
clc;

addpath(['.' filesep 'functions']);
addpath(['.' filesep 'maps']);

%n_exp = input("number of experiments: "); %number of experiments to be performed
n_exp = 1; %number of experiments to be performed
N_robots = [10 ]; %number of robots for experiments

B = load_map('ht_chantry.map');
[robotPts, ~] = loadAllScens('ht_chantry-random-.scen', 25);
[robotPts2, ~] = loadAllScens('ht_chantry-even-.scen', 25);
robotPts = [robotPts; robotPts2];

% Convert coordinates from (0,0 top-left) → (0,0 bottom-left)
[H, W] = size(B);

for i = 1:numel(robotPts)
    p = robotPts{i};
    % y_cartesian = (height - 1) - y_original
    p(:,2) = (H - 1) - p(:,2);
    robotPts{i} = p;
end

% 1) Convertir todo de golpe a una matriz 2×(2N)
M = cell2mat(robotPts.');   % transponer a 1×N para concatenar horizontalmente

% 2) Separar inicios y finales (cada par de columnas es [x y])
starts = [ M(1,1:2:end).',  M(1,2:2:end).' ];   % N×2  [xS yS]
goals  = [ M(2,1:2:end).',  M(2,2:2:end).' ];   % N×2  [xG yG]

% (opcional) Forzar enteros si procede:
% starts = round(starts); goals = round(goals);

% 3) Quitar filas con NaN/Inf (por seguridad)
starts = starts(all(isfinite(starts),2), :);
goals  = goals (all(isfinite(goals ),2), :);

% 4) Eliminar duplicados manteniendo el orden de aparición
[starts_u, iaS] = unique(starts, 'rows', 'stable');
[goals_u,  iaG] = unique(goals,  'rows', 'stable');

% 5) Volver a cell {K×1} con cada fila [x y]
initPts = mat2cell(starts_u, ones(size(starts_u,1),1), 2);
goalPts = mat2cell(goals_u,  ones(size(goals_u,1),1), 2);

fprintf('Generated %d unique initial points and %d unique final points.\n', ...
        numel(initPts), numel(goalPts));

% Flip B vertically to match Cartesian coordinates
B = flipud(B);

T = build_topology_from_B(B);                     % T.adj (NxN, sparse), T.free, T.size, ...
T.map2D = B;
% Extract indices of free cells
fullAdj = T.adj;
N       = size(fullAdj, 1);
if isfield(T, 'rem_cells') && ~isempty(T.rem_cells)
    freeIdx = T.rem_cells(:);
else
    error('Invalid T structure: T.rem_cells is missing or empty.');
end
% Compute (x,y) coordinates of the center of each free cell
for k = 1:numel(freeIdx)
    [r, c] = ind2sub([H, W], freeIdx(k));
    T.centr{k} = [c - 0.5, r - 0.5];   % <-- center, not corner
end

% Subgraph on free cells (sparse, binary, symmetric, zero-diagonal)
adj = spones(fullAdj(freeIdx, freeIdx));         % ensure 0/1
adj = adj | adj.';                               % force symmetry (undirected)
adj = adj - diag(diag(adj));                     % clear diagonal

%% --- Petri net from adjacency on free cells ---
[Pre, Post] = construct_PN(adj);
[nplaces, ntrans] = size(Pre);

fprintf(1,"\nThe Petri net has %i places and %i transitions.\n",nplaces,ntrans);

% Index mappings
invMap          = zeros(N, 1);                   % original -> reduced (0 if obstacle)
invMap(freeIdx) = 1:numel(freeIdx);
fwdMap          = freeIdx;                       % reduced -> original (linear index in B)


plot_animation = 0;%input("Do you want to plot the environment and the trajectories? (1 - yes, 0 - no)\n");

flag_ILP = 1;%input("Do you want to solve also the ILP formulation? This might take a while... (1 - yes, 0 - no)\n");

%%
for i = 1 : numel(N_robots)
    N_r = N_robots(i); % Set the number of robots for the current iteration
    if (i > 1)
        fprintf(1,'Solved %i from %i experiments.',success,n_exp);
    end
    success = 0; % Initialize success counter for the current number of robots
    for exp=1:n_exp
        fprintf(1,"\n\n=======================================\n");
        fprintf(1,"Experiment number %i (%i robots)\n",exp,N_r);
        fprintf(1,"=======================================\n");
        %% --- Randomly select N_r start/goal pairs ---
        nTotal = min(numel(initPts),numel(goalPts));
        if N_r > nTotal
            warning('Number of robots (%d) exceeds the available pairs (%d). Using all.', N_r, nTotal);
            N_r = nTotal;
        end
        rng('shuffle');                                  % different sample each run
        selIdxStart      = randperm(numel(initPts), N_r);
        selectedStart = initPts(selIdxStart);
        selIdxFin      = randperm(numel(goalPts), N_r);
        selectedFin = goalPts(selIdxFin);

        fprintf('Selected %d start and goal points.\n', N_r);

        [m0, mf, idxStart, idxGoal] = initial_marking_multi_new(selectedStart,selectedFin, B, invMap, nplaces);

        T.props = idxGoal;

        if plot_animation
           %plot_environment_new(selectedPts, T.map2D, T);
           plot_environment_new_SG(selectedStart, selectedFin, T.map2D, T);
        end

        [optVal, flag] = solve_LPs_collision_avoidance_CM(Post,Pre,mf,m0,T,flag_ILP,plot_animation);
        if flag, success = success + 1; end

        sim(exp).optim = optVal;
        sim(exp).flag  = flag;
        sim(exp).m0    = m0;
        sim(exp).mf    = mf;
        sim(exp).T     = T;
        sim(exp).success = flag;
    end
    save(sprintf('simulations_TAPF_%drobots.mat', N_r), 'sim', '-v7.3');
    fprintf(1,'\n');
    clear sim;
end

