clear;
close all;
clc;

addpath(['.' filesep 'functions']);

%% input data
fprintf(2,"==========================================================================================");
fprintf(2,"\nEfficient path planning method for task allocation algorithm for Boolean specifications\n");
fprintf(2,"==========================================================================================\n");
fprintf("The initial positions of the robots and the environment are randomly generated.\n");

flag_SP = 0;
case_t  = -1;
fprintf("Please select the type of case solved\n");
while case_t~=1 && case_t~=2
    case_t = input("reachability case (1) or Boolean specification (2): ");
    if case_t==1
        fprintf(2,"Reachability case option\n");
        fprintf("\t - robots are represented with red;\n");
        fprintf("\t - regions are represented with blue (if intermediary markings are used for collision avoidance,\n");
        fprintf("\t\t then the already visited regions are colored with green);\n");
        fprintf("\t - obstacles are represented with gray;\n");
    elseif case_t==2
        fprintf(2,"Boolean specification option - the Boolean formula in CNF form is randomly generated\n");
        fprintf("\t - robots are represented with red;\n");
        fprintf("\t - regions are represented with blue;\n");
        fprintf("\t - obstacles are represented with gray;\n");
    end
end

N_r = input("Number of robots for the experiment: ");
N_p = inf;
while N_p>N_r
    if case_t==1
        N_p = input("Number of tasks for the experiment (it should be less or equal with the number of robots): ");
    else
        N_p = input("Number of tasks for the experiment: ");
        break;
    end
end

flag_ILP = -1;
while flag_ILP~=0 && flag_ILP~=1
    flag_ILP = input("Do you want to compute the solution of the ILP problem? (1 - yes, 0 - no): ");
end

N_exp = 1; %number of experiments

env_bounds = [0,200,0,200]; %the environment bounds
env_limit = 200;
obs_size = 10; %size for the grid
n_cells = ((env_bounds(2)-env_bounds(1))/obs_size) * ((env_bounds(4)-env_bounds(3))/obs_size); %number of cells of the grid
[C,adj] = grid_decomposition_regions(env_bounds,obs_size);
T.Vert = C;
T.adj = adj;
for i = 1 : length(T.Vert)
    T.centr{i} = mean(T.Vert{i},2)';
end

for exp=1:N_exp
    %% generate new environment (at every 10 experiments)
    % number of obstacles increases with each generated environment (based on the experiment number)
    % initial positions of the robots are in the maximum connected component
    % regions of interest are in the maximum connected component
    fprintf(1,'====================================\n');
    fprintf(1,'     Experiment number %i\n',exp);
    fprintf(1,'====================================\n');

    if mod(exp-1,10)==0
        N_o = randi([round(0.2*(0.5+0.0005*exp)*n_cells) round(0.3*(0.5+0.0005*exp)*n_cells)]); %number of obstacles
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
    end

    if case_t==1
        flag = 0;
        while flag == 0 %repeat until the first problem is feasible
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

            if N_r==N_p
                mf = zeros(size(Pre,1),1);
            end
            while (length(final_robots) < N_p)
                rob = randi([round(0.5*size(Pre,1))+1 size(Pre,1)],1,N_p-length(final_robots));
                final_robots = unique([final_robots rob],'stable');
            end
            for i = 1 : N_p
                T.props(i) = T.rem_cell(final_robots(i));
                V(i,final_robots(i)) = 1;
                if N_r==N_p
                    mf(final_robots(i)) = mf(final_robots(i)) + 1;
                end
            end
            % save the environment data in the simulation structure
            if N_r~=N_p
                [simulation(exp).optim,flag] = solve_LPs_TR(Post,Pre,m0,V,T,env_limit,1,flag_ILP);
            else
                [simulation(exp).optim,flag] = solve_LPs(Post,Pre,m0,mf,T,env_limit,1,flag_ILP);
            end
        end %repeat until the first problem is feasible
    else
        N_c = N_r;
        [At,bt] = generateBoolean_formula3(N_c,N_p,N_r);

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
            [simulation(exp).optim,flagRX] = solve_LPs_boolean4_randomX(At,bt,Post,Pre,V,m0,T,env_limit,1,flag_ILP,flag_SP);

        end %repeat until the first problem is feasable
    end

    %save simulation data
    simulation(exp).n_R = N_r;
    simulation(exp).n_P = N_p;
    simulation(exp).env.n_o = N_o;
    simulation(exp).Post = Post;
    simulation(exp).Pre = Pre;
    simulation(exp).m0 = m0;
    simulation(exp).T = T;
    simulation(exp).env_limit = env_limit;

end
