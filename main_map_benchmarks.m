close all;
clear;
clc;

addpath(['.' filesep 'functions']);

%% parameters
fprintf(1,"==========================================================================================");
fprintf(1,"\nEfficient path planning method for task allocation algorithm for Boolean specifications\n");
fprintf(1,"==========================================================================================\n");
fprintf(1,"main program for benchmarks from [1] -  number of tasks is equal with the number of robots.\n\n");
fprintf("[1] Multi-Agent Pathfinding: Definitions, Variants, and Benchmarks\n");
fprintf("Roni Stern and Nathan R. Sturtevant and Ariel Felner and Sven Koenig and Hang Ma and Thayne T. Walker and \n");
fprintf("\t Jiaoyang Li and Dor Atzmon and Liron Cohen and T. K. Satish Kumar and Eli Boyarski and Roman Bartak\n");
fprintf("Symposium on Combinatorial Search (SoCS), 2019, 151-158\n\n");

fprintf("\t - initial positions of the robots are represented with red triangles;\n");
fprintf("\t - regions are represented with yellow diamonds;\n");
fprintf("\t - obstacles are represented with black;\n");

obs_size = 1;
% n_exp = 1; %number of experiments to be performed
n_exp = input("number of experiments: "); %number of experiments to be performed
% N_r = 100; %number of robots
N_r = input("number of robots: "); %number of robots

plot_animation = input("Do you want to plot the environment and the trajectories? (1 - yes, 0 - no)\n");

flag_ILP = input("Do you want to solve also the ILP formulation? This might take a while... (1 - yes, 0 - no)\n");

scenario = input("Please choose the scenario (1 - room-32-32-4, 2 - random-32-32-20, 3 - den312d, 4 - ht_chantry) ");

if scenario==1
    map_file = 'room-32-32-4.map';
    scenario_file = 'room-32-32-4-random-1-25.xlsx';
    % R = [10 50 100 200 300];
elseif scenario==2
    map_file = 'random-32-32-20.map';
    scenario_file = 'random-32-32-20-random-1-25.xlsx';
    % R = [10 50 100 200 300 400];
elseif scenario==3
    map_file = 'den312d.map';
    scenario_file = 'den312d-random-1-25.xlsx';
    % R = [10 50 100 200 300 400 500 600 700 800 900 1000 1100 1200];
elseif scenario==4
    map_file = 'ht_chantry.map';
    scenario_file = 'ht_chantry-random-1-25.xlsx';
    % R = [10 50 100 200 300 400 500 600 700 800 900 1000 1100 1500 2000];
end

%% 
T = grid_decomposition_regions(map_file, obs_size);
env_limit = size(T.map2D);

fullAdj = T.adj;
freeMask = true(size(fullAdj,1),1);
freeMask(T.obstacles) = false;
adj = fullAdj(freeMask,freeMask);

[Pre,Post] = construct_PN(adj);

invMap = zeros(numel(freeMask),1);
invMap(T.rem_cell) = 1:numel(T.rem_cell);

scen        = readmatrix(scenario_file);
sx          = scen(:,5); sy = scen(:,6);
gx          = scen(:,7); gy = scen(:,8);

linStartAll = sub2ind(env_limit, sy+1, sx+1);
linGoalAll  = sub2ind(env_limit, gy+1, gx+1);

okStart     = T.map2D(linStartAll);
okGoal      = T.map2D(linGoalAll);
maskBoth    = okStart & okGoal;

validStarts = unique(linStartAll(maskBoth));
validGoals  = unique(linGoalAll(maskBoth));

% if numel(validStarts) < max(R) || numel(validGoals) < max(R)
%     error('Insuficientes celdas únicas libres: starts=%d, goals=%d, maxNeeded=%d', ...
%         numel(validStarts), numel(validGoals), max(R));
% end

success = 0;
for exp=1:n_exp
    fprintf("\nExperiment number %i (%i robots)",exp,N_r);

    flag = 0;
    while flag==0
        selS = randperm(numel(validStarts), N_r);
        goals = setdiff(validGoals,validStarts(selS)); %startSamples != goalSamples
        selG = randperm(numel(goals),  N_r);
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

        if exp==1 && plot_animation
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

% if scenario==1
%     save(sprintf('room-32-32-4_%drobots.mat', N_r), 'sim', '-v7.3');
% elseif scenario==2
%     save(sprintf('random-32-32-20_%drobots.mat', N_r), 'sim', '-v7.3');
% elseif scenario==3
%     save(sprintf('den312d_%drobots.mat', N_r), 'sim', '-v7.3');
% elseif scenario==4
%     save(sprintf('ht-chantry_%drobots.mat', N_r), 'sim', '-v7.3');
% end