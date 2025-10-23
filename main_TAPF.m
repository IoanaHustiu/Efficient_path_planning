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

obs_size = 1;
n_exp = 20; %number of experiments to be performed
%n_exp = input("number of experiments: "); %number of experiments to be performed
% N_r = 100; %number of robots
N_robots = [100 250 500 750 1000 1500]; %number of robots

B = load_map('Paris_1_256.map');
robotPts = loadParis256AllScens(".\maps");
%B = load_map('lgt601d.map');
%robotPts = loadSingleScenFile('lgt601d_map.scen');
% Convert coordinates from (0,0 top-left) → (0,0 bottom-left)
[H, W] = size(B);

for i = 1:numel(robotPts)
    p = robotPts{i};
    % y_cartesian = (height - 1) - y_original
    p(:,2) = (H - 1) - p(:,2);
    robotPts{i} = p;
end
% Flip B vertically to match Cartesian coordinates
B = flipud(B);

T = build_topology_from_B(B);                     % T.adj (NxN, sparse), T.free, T.size, ...
T.map2D = B;
% Extract indices of free cells
fullAdj = T.adj;
N       = size(fullAdj, 1);
if isfield(T, 'free') && ~isempty(T.free)
    freeIdx = T.free(:);
else
    error('Invalid T structure: T.free is missing or empty.');
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
        fprintf(1,"\n=======================================\n");
        fprintf(1,"Experiment number %i (%i robots)\n",exp,N_r);
        fprintf(1,"=======================================\n");
        %% --- Randomly select N_r start/goal pairs ---
        nTotal = numel(robotPts);
        if N_r > nTotal
            warning('Number of robots (%d) exceeds the available pairs (%d). Using all.', N_r, nTotal);
            N_r = nTotal;
        end
        rng('shuffle');                                  % different sample each run
        selIdx      = randperm(nTotal, N_r);
        selectedPts = robotPts(selIdx);

        fprintf('Selected %d out of %d start/goal pairs.\n', N_r, nTotal);

        [m0, mf, idxStart, idxGoal] = initial_marking_multi(selectedPts, B, invMap, nplaces);

        T.props = idxGoal;

        %if plot_animation
        %   plot_environment(startSamples, goalSamples, T.map2D, env_limit, T.Vert);
        %   title(sprintf('Start/Goal for %d robots', N_r));
        %end

        [optVal, flag] = solve_LPs_collision_avoidance(Post,Pre,mf,m0,T,flag_ILP,plot_animation);
        if flag, success = success + 1; end

        sim(exp).optim = optVal;
        sim(exp).flag  = flag;
        sim(exp).m0    = m0;
        sim(exp).mf    = mf;
        sim(exp).T     = T;
        sim(exp).success = flag;
    end
    save(sprintf('simulations_TAPF_%drobots.mat', N_r), 'sim', '-v7.3');
    clear sim;
end


% save(sprintf('simulations_reachability_%drobots.mat', N_r), 'sim', '-v7.3');
% clear sim;