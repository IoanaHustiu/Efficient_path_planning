close all;
clear;
clc;

addpath(['.' filesep 'functions']);

%% parameters
fprintf(1,"==========================================================================================");
fprintf(1,"\nEfficient path planning method for task allocation algorithm for Boolean specifications\n");
fprintf(1,"==========================================================================================\n");
fprintf(1,"Main program for smart plant scenario\n");
fprintf("The initial positions of the robots and the environment are randomly generated.\n");

fprintf("\t - initial positions of the robots are represented with red triangles;\n");
fprintf("\t - regions are represented with yellow diamonds;\n");
fprintf("\t - obstacles are represented with black;\n");

obs_size = 1;
N_r = 20; %number of robots
% N_r = input("number of robots: "); %number of robots
map = [N_r N_r];

plot_animation = input("Do you want to plot the environment and the trajectories? (1 - yes, 0 - no)\n");

flag_ILP = input("Do you want to solve also the ILP formulation? This might take a while... (1 - yes, 0 - no)\n");

nj1 = 2; %jobs type1
nj2 = 2; %jobs type2
nj3 = 2; %jobs type3
nj4 = 4; %jobs type4
N_p = nj1 + nj2 + nj3 + nj4; %manufacturing plant

fprintf("Number of robots: %i\n",N_r);
fprintf("\t - the robots that are required to move are represented with red;\n");
fprintf("\t - the robots that are not required to move are represented with green;\n");
fprintf("Number of tasks: %i, from which\n",N_p);
fprintf("\t\t - type 1 jobs: %i (should be 100%% fulfilled) - the blue cells\n",nj1);
fprintf("\t\t - type 2 jobs: %i (should be 100%% fulfilled) - the cyan cells\n",nj2);
fprintf("\t\t - type 3 jobs: %i (should be at least 50%% fulfilled) - the dark green cells\n",nj3);
fprintf("\t\t - type 4 jobs: %i (should be at least 75%% fulfilled) - the yellow cells\n",nj4);

%%
map_size = [map(1),map(2),round(0.25*map(1)*map(2))];

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

% if numel(validStarts) < N_r || numel(validGoals) < N_r
%     error('Insuficientes celdas únicas libres: starts=%d, goals=%d, maxNeeded=%d', ...
%           numel(validStarts), numel(validGoals), N_r);
% end

%Generate Boolean goal for smart plant scenario
bloque1 = -eye(nj1); %job type 1 - 100%;
bloque2 = -eye(nj2); %job type 2 - 100%;
c_job3 = nchoosek(1:nj3,6); %job type 3 - 50%;
bloque3 = zeros(size(c_job3,1),nj3);
for i=1:size(c_job3,1)
    bloque3(i,c_job3(i,:)) = -1;
end
c_job4 = nchoosek(1:nj4,4); %job type4 - 75%;
bloque4 = zeros(size(c_job4,1),nj4);
for i=1:size(c_job4,1)
    bloque4(i,c_job4(i,:)) = -1;
end
At = [bloque1 zeros(size(bloque1,1),size(bloque2,2)) zeros(size(bloque1,1),size(bloque3,2)) zeros(size(bloque1,1),size(bloque4,2));...
    zeros(size(bloque2,1),size(bloque1,2)) bloque2 zeros(size(bloque2,1),size(bloque3,2)) zeros(size(bloque2,1),size(bloque4,2));...
    zeros(size(bloque3,1),size(bloque1,2)) zeros(size(bloque3,1),size(bloque2,2)) bloque3 zeros(size(bloque3,1),size(bloque4,2));...
    zeros(size(bloque4,1),size(bloque1,2)) zeros(size(bloque4,1),size(bloque2,2)) zeros(size(bloque4,1),size(bloque3,2)) bloque4];
bt = (sum(At'==1)-1)';

success = 0;
flag_ILP = 1;
plot_animation = 1;


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
% mf = zeros(size(Pre,1),1);
% mf(goalIdx)  = 1;

T.props = goalIdx;
T.nj = [nj1 nj2 nj3 nj4];

if plot_animation
    plot_environment_smartPlant(T.nj, startSamples, goalSamples, T.obstacles, env_limit, T.Vert);
    title(sprintf('Start/Goal for %d robots', N_r));
end

[optVal, flag] = solve_LPs_collision_avoidance_boolean(Post,Pre,At,bt,m0,T,flag_ILP,plot_animation);
if flag, success = success + 1; end

sim.optim = optVal;
sim.flag  = flag;
sim.m0    = m0;
sim.T     = T;
sim.success = flag;
sim.goalSamples = goalSamples;
sim.goals = goalIdx;



% success_rate = success / n_exp * 100;
% fprintf('\nSuccess rate for %d robots (smart manufacturing plant scenario): %.2f%%\n', N_r, success_rate);

% save(sprintf('simulations_%drobots_smart_plant.mat', N_r), 'sim', '-v7.3');
% clear sim;