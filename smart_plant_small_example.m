clear;
close all;
warning on;
clc;

addpath(['.' filesep 'functions']);

%% input data for the simulation
fprintf(2,"==========================================================================================");
fprintf(2,"\nEfficient path planning method for task allocation algorithm for Boolean specifications\n");
fprintf(2,"==========================================================================================\n");

fprintf("SMART PLANT EXAMPLE\n");
flag_ILP = input("Do you want to compute the solution of the ILP problem? It might take a while... (1 - yes, 0 - no): ");

N_exp = 1; %number of experiments
flag_SP = 1;

env_bounds = [0,300,0,300]; %the environment bounds
env_limit = 300;
obs_size = 10; %size for the grid
N_r = 50;
nj1 = 10; %jobs type1
nj2 = 15; %jobs type2
nj3 = 10; %jobs type3
nj4 = 15; %jobs type4
N_p = nj1 + nj2 + nj3 + nj4; %manufacturing plant

fprintf("Number of robots: %i\n",N_r);
fprintf("Number of tasks: %i, from which\n",N_p);
fprintf("\t\t - type 1 jobs: %i (should be 100%% fulfilled)\n",nj1);
fprintf("\t\t - type 2 jobs: %i (should be 100%% fulfilled)\n",nj2);
fprintf("\t\t - type 3 jobs: %i (should be at least 50%% fulfilled)\n",nj3);
fprintf("\t\t - type 4 jobs: %i (should be at least 75%% fulfilled)\n",nj4);

n_cells = ((env_bounds(2)-env_bounds(1))/obs_size) * ((env_bounds(4)-env_bounds(3))/obs_size); %number of cells of the grid
[C,adj] = grid_decomposition_regions(env_bounds,obs_size);
T.Vert = C;
T.adj = adj;
T.nj = [nj1 nj2 nj3 nj4];
for i = 1 : length(T.Vert)
    T.centr{i} = mean(T.Vert{i},2)';
end

for exp=1:N_exp
    %% generate new environment (at every 10 experiments)
    % number of obstacles increases with each generated environment (based on the experiment number)
    % initial positions of the robots are in the maximum connected component
    % regions of interest are in the maximum connected component
    fprintf(1,'\n====================================\n');
    fprintf(1,'     Experiment number %i\n',exp);
    fprintf(1,'====================================\n');

    if mod(exp-1,10)==0
        N_o = randi([round(0.2*(0.5+0.0005*exp)*n_cells) round(0.3*(0.5+0.0005*exp)*n_cells)]); %number of obstacles
        % N_o = 0;
        obstacles = sort(random_env_generator(N_o,n_cells));
        adj = T.adj;
        for i = length(obstacles):-1:1
            adj(obstacles(i),:) = [];
            adj(:,obstacles(i)) = [];
        end
        [Pre,Post] = construct_PN(adj);
        C = Post-Pre;
        T.obstacles = obstacles;
        T.rem_cell = setdiff(1:n_cells,obstacles);

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

    end

    flagRX = 0;
    while flagRX == 0 %repeat until the first problem is feasible
        initial_robots = [];
            m0 = zeros(size(Pre,1),1);
            while (length(initial_robots) < N_r)
                rob = randi([1 round(0.5*size(Pre,1))],1,N_r-length(initial_robots));
                initial_robots = unique([initial_robots rob],'stable');
            end
            for i = 1 : length(initial_robots)
                m0(initial_robots(i)) = m0(initial_robots(i)) + 1;
            end

            final_robots = [];
            T.props = zeros(1,N_p);
            V = zeros(N_p,size(C,1));

            while (length(final_robots) < N_p)
                rob = randi([round(0.5*size(Pre,1))+1 size(Pre,1)],1,N_p-length(final_robots));
                final_robots = unique([final_robots rob],'stable');
            end
            for i = 1 : N_p
                T.props(i) = T.rem_cell(final_robots(i));
                V(i,final_robots(i)) = 1;
            end
            [simulation(exp).optimRX,flagRX] = solve_LPs_boolean4_randomX(At,bt,Post,Pre,V,m0,T,env_limit,1,flag_ILP,flag_SP);

    end %repeat until the first problem is feasable

    simulation(exp).n_R = N_r;
    simulation(exp).n_P = N_p;
    simulation(exp).env.n_o = N_o;
    simulation(exp).Post = Post;
    simulation(exp).Pre = Pre;
    simulation(exp).m0 = m0;
    simulation(exp).T = T;
    simulation(exp).env_limit = env_limit;
end
