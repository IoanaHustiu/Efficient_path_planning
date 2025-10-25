close all;
clear;
clc;

addpath(['.' filesep 'functions']);
addpath(['.' filesep 'maps']);

%% parameters
fprintf(1,"==========================================================================================");
fprintf(1,"\nTask Allocation and Path Findings algorithm\n");
fprintf(1,"==========================================================================================\n");
fprintf(1,"main program for TAPF -  number of tasks is equal with the number of robots.\n");
fprintf("The initial positions of the robots and the environment are randomly generated from .scen file.\n");

fprintf("\t - initial positions of the robots are represented with red triangles;\n");
fprintf("\t - regions are represented with yellow diamonds;\n");
fprintf("\t - obstacles are represented with black;\n");

%n_exp = input("number of experiments: "); %number of experiments to be performed
n_exp = 10; %number of experiments to be performed
%N_robots = [500 20 30 40 50 100 250 500 750 1000 1250 1500 1750 2000 2500 2750 3000 3250 3500 3750 4000]; %number of robots for experiments
N_robots = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16]; %number of robots for experiments

B = load_map('warehouse-10-20-10-2-4.map');
[robotPts, scenTbl] = loadAllScens('warehouse-10-20-10-2-1-random-.scen', 25);

% Convert coordinates from (0,0 top-left) → (0,0 bottom-left)
[H, W] = size(B);

for i = 1:numel(robotPts)
    p = robotPts{i};
    % y_cartesian = (height - 1) - y_original
    p(:,2) = (H - 1) - p(:,2);
    robotPts{i} = p;
end

% --- Extrae todos los puntos de inicio y final
starts = cellfun(@(p) p(1,:), robotPts, 'UniformOutput', false);
goals  = cellfun(@(p) p(2,:), robotPts, 'UniformOutput', false);

starts = vertcat(starts{:});
goals  = vertcat(goals{:});

% --- Selección según coordenada x
mask_init = (starts(:,1) < 25);
mask_goal = (goals(:,1)  > 28) & (goals(:,1)  < 130);

% --- Puntos iniciales y finales candidatos
initPts_all = starts(mask_init, :);
goalPts_all = goals(mask_goal, :);

% --- Asegurar que no hay puntos repetidos
[~, ia] = unique(initPts_all, 'rows');
initPts = num2cell(initPts_all(ia, :), 2);

[~, ig] = unique(goalPts_all, 'rows');
goalPts = num2cell(goalPts_all(ig, :), 2);

% --- Eliminar solapamiento entre conjuntos (por si acaso)
initMat = vertcat(initPts{:});
goalMat = vertcat(goalPts{:});

[~, ia] = setdiff(initMat, goalMat, 'rows', 'stable');
[~, ig] = setdiff(goalMat, initMat, 'rows', 'stable');

initPts = initPts(ia);
goalPts = goalPts(ig);

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

%plot_environment_new([], T.map2D, T);

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
           plot_environment_new_SG(selectedStart, selectedFin, T.map2D, T);
        end

        [optVal, flag] = solve_LPs_collision_avoidance(Post,Pre,mf,m0,T,flag_ILP,plot_animation);
        if flag, success = success + 1; end

        sim(exp).optim = optVal;
        sim(exp).flag  = flag;
        sim(exp).m0    = m0;
        sim(exp).mf    = mf;
        sim(exp).T     = T;
        sim(exp).success = flag;
    end
    save(sprintf('simulations_TAPF_aisle_%drobots.mat', N_r), 'sim', '-v7.3');
    clear sim;
end
