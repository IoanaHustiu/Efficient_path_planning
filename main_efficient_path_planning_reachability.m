clear;
close all;
warning on;
clc;

addpath(['.' filesep 'functions']);

%% input data for the simulation
N_exp = 10; %number of experiments
N_r = 10; %number of robots
N_p = 10; %numbers of regions of interest
flagBoolean = 1; %flag=1 -> Boolean formula, flag = 0 -> reachability

reg_edges={'r','b','g','c','m','k','y'};    %colors of proposition boundaries
rob_color={'r','b','g','c','m','k','y'};    %colors of robots

env_bounds = [0,100,0,100]; %the environment bounds
env_limit = 100;
obs_size = 10; %size for the grid
n_cells = ((env_bounds(2)-env_bounds(1))/obs_size) * ((env_bounds(4)-env_bounds(3))/obs_size); %number of cells of the grid

for exp=1:N_exp
    %% generate new environment (at every 10 experiments)
    % number of obstacles increases with each generated environment (based on the experiment number)
    % initial positions of the robots are in the maximum connected component
    % regions of interest are in the maximum connected component
    fprintf(1,'Experiment number %i\n',exp);

    if mod(exp-1,10)==0
        N_o = randi([round(0.2*(0.5+0.005*exp)*n_cells) round(0.3*(0.5+0.005*exp)*n_cells)]); %number of obstacles
        [regions,obstacles,o_cells,initial_point,~,adj,rem_cells] = random_env_generator(N_p,N_o,obs_size,N_r,env_bounds);

        if flagBoolean
            N_c = randi([1,2^N_p],1);
            [At,bt] = generateBoolean_formula(N_c,N_p);
%         At = [-1 0 0 0 0;0 0 -1 0 0;0 0 0 0 -1;0 0 0 -1 0];
%         N_c = size(At,1);
%         bt = (sum(At'==1)-1)';
        end

    else %keep the position of the obstacles and of the regions of interest and change the initial position of the robots
        random = [];
        x = (env_bounds(2)-env_bounds(1))/obs_size;
        y = (env_bounds(4)-env_bounds(3))/obs_size;
        while length(random)<N_r
            rob = randi([1 round(0.5*length(rem_cells))],1,N_r-length(random));
            random = unique([random rem_cells(rob)],'stable');
        end

        for i=1:N_r
            cells = random(i);
            pos_x = mod(cells-1,x);
            pos_y = floor((cells-1)/x);
            initial_point{i}=mean(obs_size*[pos_x pos_x+1 pos_x+1 pos_x; pos_y pos_y pos_y+1 pos_y+1],2);
        end
    end

    %% create environment structure
    T = create_partition_obstacles(regions,o_cells,env_bounds,initial_point);    %partition, adjacency, observations, initial cells of robots
    %% create RMPN model
    [Pre,Post,~] = construct_PN(T);
    C = Post-Pre;

    V = zeros(N_p,size(C,1));
    for i=1:N_p
        V(i,T.props{i}) = 1;
    end

    m0 = zeros(size(Pre,1),1);
    for i = 1 : length(T.R0)
        m0(T.R0(i)) = m0(T.R0(i)) + 1;
    end

    if flagBoolean==0
        mf = zeros(size(Pre,1),1);
        for i = 1 : length(T.props)
            mf(T.props{i}) = mf(T.props{i}) + 1;
        end
    end

    % save the environment data in the simulation structure
    simulation(exp).n_R = N_r;
    simulation(exp).n_P = N_p;
    simulation(exp).env.n_o = N_o;
    simulation(exp).Post = Post;
    simulation(exp).Pre = Pre;
    simulation(exp).m0 = m0;
    simulation(exp).T = T;
    simulation(exp).env_limit = env_limit;

    if flagBoolean
        simulation(exp).optim = solve_LPs_boolean(At,bt,Post,Pre,V,m0,T,env_limit,0);
    else
        simulation(exp).mf = mf;
        simulation(exp).optim = solve_LPs(Post,Pre,m0,mf,T,env_limit,0);
    end

end
