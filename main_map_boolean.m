close all;
clear;
clc;

addpath(['.' filesep 'functions']);

%% parameters
fprintf(1,"==========================================================================================");
fprintf(1,"\nEfficient path planning method for task allocation algorithm for Boolean specifications\n");
fprintf(1,"==========================================================================================\n");
fprintf(1,"main program for boolean-based goal\n");
fprintf("The initial positions of the robots, the environment and the Boolean-based goal are randomly generated.\n");

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
        map_size = [map(1),map(2),round(0.25*map(1)*map(2))];

        %
        T = grid_decomposition_regions_environment(map_size, obs_size);
        env_limit  = size(T.map2D,1);

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

        N_p = randi([round(0.75*N_r),round(1.25*N_r)]);
        N_c = N_p;

        %Generate Boolean goal
        At = -true(N_c,N_p);  % Start with all -1
        % Randomly choose indices to set to 0
        k = randi([1,N_c * N_p]);
        idx = randperm(N_c * N_p, k);
        At(idx) = 0;

        At = unique(At,'rows','stable'); %eliminate duplicate rows
        At = At(~all(At==0,2),:); %eliminate rows containing only entries equal with 0
        N_c = size(At,1); %update the number of conjunctions (if necessary)
        cols_all_zero = all(At == 0, 1);  % Logical row vector: 1 if entire column is zeros
        zero_col_ind = find(cols_all_zero);
        At(randi([1,N_c],1,length(zero_col_ind)),zero_col_ind) = -1;

        bt = (sum(At' == 1)-1)';
    end

    success = 0;
    flag = 0;

    while flag==0
        selS = randperm(numel(validStarts), N_r);
        validGoals = setdiff(validGoals,validStarts(selS)); %startSamples != goalSamples
        selG = randperm(numel(validGoals),  N_p);
        startSamples = validStarts(selS);
        goalSamples  = validGoals(selG);

        startIdx = invMap(startSamples);
        goalIdx  = invMap(goalSamples);

        assert(numel(unique(startIdx))==N_r, 'Duplicado en startIdx');
        assert(numel(unique(goalIdx ))==N_p, 'Duplicado en goalIdx');

        m0 = zeros(size(Pre,1),1);
        m0(startIdx) = 1;

        T.props = goalIdx;

        if exp==1 && plot_animation
            plot_environment(startSamples, goalSamples, T.map2D, env_limit, T.Vert);
            title(sprintf('Start/Goal for %d robots', N_r));
        end

        [optVal, flag] = solve_LPs_collision_avoidance_boolean(Post,Pre,At,bt,m0,T,flag_ILP,plot_animation);
        if flag, success = success + 1; end
    end

    sim(exp).optim = optVal;
    sim(exp).flag  = flag;
    sim(exp).m0    = m0;
    sim(exp).T     = T;
    sim(exp).success = flag;
end

% success_rate = success / n_exp * 100;
% fprintf('\nSuccess rate for %d robots: %.2f%%\n', N_r, success_rate);

% save(sprintf('simulations_boolean_%drobots.mat', N_r), 'sim', '-v7.3');
% clear sim;