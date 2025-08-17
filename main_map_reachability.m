close all;
clear;
clc;

addpath(['.' filesep 'functions']);

%% parameters
fprintf(1,"==========================================================================================");
fprintf(1,"\nEfficient path planning method for task allocation algorithm for Boolean specifications\n");
fprintf(1,"==========================================================================================\n");
fprintf(1,"main program for reachability case -  number of tasks is equal with the number of robots.\n");
fprintf("The initial positions of the robots and the environment are randomly generated.\n");

fprintf("\t - initial positions of the robots are represented with red triangles;\n");
fprintf("\t - regions are represented with yellow diamonds;\n");
fprintf("\t - obstacles are represented with black;\n");

obs_size = 1;
% n_exp = 1; %number of experiments to be performed
n_exp = input("number of experiments: "); %number of experiments to be performed
% N_r = 100; %number of robots
N_r = input("number of robots: "); %number of robots
map = [N_r N_r];

plot_animation = input("Do you want to plot the environment and the trajectories? (1 - yes, 0 - no)\n");

flag_ILP = input("Do you want to solve also the ILP formulation? This might take a while... (1 - yes, 0 - no)\n");

%%
for exp=1:n_exp
    fprintf(2,"\nExperiment number %i (%i robots)\n",exp,N_r);

    if mod(exp-1,10)==0
        map_size      = [map(1),map(2),round(0.25*map(1)*map(2))];

        %
        T = grid_decomposition_regions_environment(map_size, obs_size);
        env_limit  = size(T.map2D);

        fullAdj = T.adj;
        freeMask = true(size(fullAdj,1),1);
        freeMask(T.obstacles) = false;
        adj = fullAdj(freeMask,freeMask);

        [Pre,Post] = construct_PN(adj);

        invMap = zeros(numel(freeMask),1);
        invMap(T.rem_cell) = 1:numel(T.rem_cell);

        %
        % Find connected components (4-connectivity)
        CC = bwconncomp(full(T.map2D), 4);

        % Get sizes of each component
        component_sizes = cellfun(@numel, CC.PixelIdxList);

        % Find index of the largest component
        [~, idx] = max(component_sizes);

        % Get the set of linear indices for that component
        largest_component_indices = CC.PixelIdxList{idx};

        linStartAll = largest_component_indices;
        linGoalAll  = largest_component_indices;

        %
        okStart     = T.map2D(linStartAll);
        okGoal      = T.map2D(linGoalAll);
        maskBoth    = okStart & okGoal;

        validStarts = unique(linStartAll(maskBoth));
        validGoals  = unique(linGoalAll(maskBoth));

        % if numel(validStarts) < max(R) || numel(validGoals) < max(R)
        %     error('Insuficientes celdas únicas libres: starts=%d, goals=%d, maxNeeded=%d', ...
        %         numel(validStarts), numel(validGoals), max(R));
        % end
    end

    success = 0;
    flag = 0;
    while flag==0
        selS = randperm(numel(validStarts), N_r);
        selG = randperm(numel(validGoals),  N_r);
        startSamples = validStarts(selS);
        goalSamples  = validGoals(selG);

        startIdx = invMap(startSamples);
        goalIdx  = invMap(goalSamples);

        assert(numel(unique(startIdx))==N_r, 'Duplicado en startIdx');
        assert(numel(unique(goalIdx ))==N_r, 'Duplicado en goalIdx');

        m0 = zeros(size(Pre,1),1);
        m0(startIdx) = 1;
        mf = zeros(size(Pre,1),1);
        mf(goalIdx)  = 1;

        T.props = goalIdx;

        if plot_animation
            plot_environment(startSamples, goalSamples, T.map2D, env_limit, T.Vert);
            title(sprintf('Start/Goal for %d robots', N_r));
        end

        [optVal, flag] = solve_LPs_collision_avoidance(Post,Pre,mf,m0,T,flag_ILP,plot_animation);
        if flag, success = success + 1; end
    end

    sim(exp).optim = optVal;
    sim(exp).flag  = flag;
    sim(exp).m0    = m0;
    sim(exp).mf    = mf;
    sim(exp).T     = T;
    sim(exp).success = flag;
end

% save(sprintf('simulations_reachability_%drobots.mat', N_r), 'sim', '-v7.3');
% clear sim;