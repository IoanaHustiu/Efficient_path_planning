clear;
close all;
warning on;
clc;

addpath(['.' filesep 'functions'],['.' filesep 'simulations']);
runSolve();

%% initialize structure for saving simulation data
% simulation.RMPN = []; %add the number of cells, obstacles, places and transitions at each experiment

% for each simulation the runtime, the cost, the number of iteration/intermediary markings and the cell capacity will be saved
simulation.completeMILP = [];
simulation.completeLP = [];
simulation.interMILP = [];
simulation.interLP = [];

%% input data for the simulation
N_exp = 2; %number of experiments
N_r = 30; %number of robots
N_p = 30; %numbers of regions of interest

simulation.env.n_R = N_r;
simulation.env.n_P = N_p;

reg_edges={'r','b','g','c','m','k','y'};    %colors of proposition boundaries
rob_color={'r','b','g','c','m','k','y'};    %colors of robots

env_bounds = [0,200,0,200]; %environment bounds
obs_size = 10; %size for the grid
n_cells = ((env_bounds(2)-env_bounds(1))/obs_size) * ((env_bounds(4)-env_bounds(3))/obs_size); %number of cells of the grid

for exp=1:N_exp
%% generate new environment (at every 10 experiments)
% number of obstacles increases with each generated environment (based on the experiment number)
% initial positions of the robots are in the maximum connected component
% regions of interest are in the maximum connected component
    fprintf('Experiment number %i\n',exp);

    if mod(exp-1,10)==0
        N_o = randi([round(0.2*(0.5+0.005*exp)*n_cells) round(0.3*(0.5+0.005*exp)*n_cells)]); %number of obstacles
%         N_o = randi([round(0.2*n_cells) round(0.3*n_cells)]); %number of obstacles
        [regions,obstacles,o_cells,initial_point,~,adj,rem_cells] = random_env_generator(N_p,N_o,obs_size,N_r,env_bounds);
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

% save the initial positions of the robots
    x0 = initial_point;

%% Boolean formula - reachability case
    At = -eye(N_p);
    bt = -ones(N_p,1);

%% create environment structure
    T = create_partition_obstacles(regions,o_cells,env_bounds,x0);    %partition, adjacency, observations, initial cells of robots

%     plot_environment_obstacles(regions,obstacles,env_bounds,reg_edges,T.Vert); %partition
%     for r=1:N_r    %plot initial positions of robots
%         plot(x0{r}(1,1),x0{r}(2,1),rob_color{mod(r-1,length(rob_color))+1},'Marker','o','LineWidth',1.5);
%     end

%% create RMPN model
    [Pre,Post,m0] = construct_PN(T);
    C = Post-Pre;

    np = size(C,1); %number of places
    nt = size(C,2); %number of transitions

% save the environment data in the simulation structure
    simulation.env.n_o(exp) = N_o;
    simulation.env.n_places(exp) = np-N_o;
    simulation.env.n_transitions(exp) = sum(sum(full(T.adj)-diag(diag(full(T.adj)))));
    simulation.env.initialPlaces(exp) = {m0};
    simulation.env.initialPositions(exp) = {x0};

% construct matrix V over the alphabet
    V = zeros(N_p,size(C,1));
    for i=1:N_p
        V(i,T.props{i}) = 1;
    end

    RobotInitPlaces = T.R0; %collect the initial places of the robots

    % compute the initial marking
%     [r0,r0Index] = groupcounts(RobotInitPlaces');
%     m0 = zeros(size(Post,1),1);
%     m0(r0Index) = r0;

    M = 1e4; %weight of the second term in the objective function

%% workspace
%     load('interMILP_unfeasible.mat');

%% (1) MILP complete problem - s and sigma as variables for reachability case (with solution (sM,sigmaM))
% x is known (has all entries equal with 1)
    x = ones(N_p,1);
% m is known (V'*x)
    mf = V'*x;

    sigma = optimvar('sigma',nt,1,'LowerBound',0,'UpperBound',inf,'Type','integer');
    s = optimvar('s',1,1,'LowerBound',0,'UpperBound',N_r,'Type','integer');
    
    prob = optimproblem;
    prob.Objective = sum(sigma) + M*s;

    prob.Constraints.cons1 = mf == m0 + C*sigma;
    prob.Constraints.cons2 = Post*sigma + m0 <= s*ones(np,1);
    
    opt = optimoptions('intlinprog','Display','none');
    intcon = 1:nt+1;

    tic;
    [sol_completeMILP,fval_completeMILP,exitflag_completeMILP,output_completeMILP] = solve(prob,'Options',opt);
    time_completeMILP = toc;

    if exitflag_completeMILP~=1 && exitflag_completeMILP~=2
        fprintf('MILP complete problem: intlinprog stopped - problem unfeasible.\n');
        return;
    else
        tol = opt.ConstraintTolerance;

        if norm(sol_completeMILP.sigma-round(sol_completeMILP.sigma))<=10*tol
            sigma_completeMILP = round(sol_completeMILP.sigma);
            s_completeMILP = round(sol_completeMILP.s);

            cost_completeMILP = sum(sigma_completeMILP);

% compute trajectories                
                rob_traj_completeMILP = compute_trajectories(T,Pre,Post,m0,sigma_completeMILP,x0);
% plot trajectories 
                plot_environment_obstacles(regions,obstacles,env_bounds,reg_edges,T.Vert); %partition
                title('relaxed LP complete problem');
                plot_trajectories(rob_traj_completeMILP,rob_color);
% plot animation                
                plot_animation(T,Pre,Post,[m0 mf],sigma_completeMILP,env_bounds);
        else
            fprintf('MILP complete problem - sigma is not integer!\n');
        end
    end

% save the data in simulation structure
    simulation.completeMILP.runtime(exp) = time_completeMILP;
    simulation.completeMILP.cost(exp) = cost_completeMILP;
    simulation.completeMILP.cellCapacity(exp) = s_completeMILP;
    simulation.completeMILP.sigma(exp) = {sigma_completeMILP};

%% (2) relaxed LP complete problem - s and sigma as variables (with solution (s*,sigma*))
% change the type of variables from integer to continuous
    sigma.Type = 'continuous';
    s.Type = 'continuous';

    opt = optimoptions('intlinprog','Display','none');
    intcon = [];

    tic;
    [sol_completeLP,fval_completeLP,exitflag_completeLP,output_completeLP] = solve(prob,'Options',opt);
    time_completeLP = toc;

    if exitflag_completeLP~=1 && exitflag_completeLP~=2
        fprintf('relaxed LP complete problem: intlinprog stopped - problem unfeasible.\n');
%         return;
    else
        cost_completeLP = Inf;
        tol = opt.ConstraintTolerance;
        unrounded_s_completeLP = sol_completeLP.s;

        if norm(sol_completeLP.s-round(sol_completeLP.s))<=10*tol
            fprintf('relaxed LP complete problem: integer cell capacity.\n');
            sigma_completeLP = round(sol_completeLP.sigma);
            s_completeLP = round(sol_completeLP.s);

            cost_completeLP = sum(sigma_completeLP);

            if cost_completeLP~=cost_completeMILP
                fprintf('Different cost between MILP and relaxed LP.\n');
            end
        else
            fprintf('relaxed LP complete problem - sigma is not integer!\n');
            
            s_bar = ceil(sol_completeLP.s);

            %solve LP2 only once since s will be integer
            sigma = optimvar('sigma',nt,1,'LowerBound',0,'UpperBound',inf,'Type','continuous');

            prob = optimproblem;
            prob.Objective = sum(sigma);

            prob.Constraints.cons1 = mf == m0 + C*sigma;
            prob.Constraints.cons2 = Post * sigma + m0 <= s_bar*ones(np,1);

            opt = optimoptions('intlinprog','Display','none');

            tic;
            [sol_completeLP2,fval_completeLP2,exitflag_completeLP2,output_completeLP2] = solve(prob,'Options',opt);
            time_completeLP = time_completeLP + toc;

            if exitflag_completeLP2~=1 && exitflag_completeLP2~=2
                fprintf('Iterative complete LP is not feasible.\n');
            else
                if norm(sol_completeLP2.sigma-round(sol_completeLP2.sigma))<=10*tol
                    fprintf('relaxed LP complete problem: integer cell capacity.\n');
                    sigma_completeLP = round(sol_completeLP2.sigma);

                    cost_completeLP = sum(sigma_completeLP);
                else
                    fprintf('relaxed LP complete problem - sigma is not integer!\n');
                end
            end

% compute trajectories                
                rob_traj_completeLP = compute_trajectories(T,Pre,Post,m0,sigma_completeLP,x0);
% plot trajectories 
                plot_environment_obstacles(regions,obstacles,env_bounds,reg_edges,T.Vert); %partition
                title('relaxed LP complete problem');
                plot_trajectories(rob_traj_completeLP,rob_color);
% plot animation                
                plot_animation(T,Pre,Post,[m0 mf],sigma_completeLP,env_bounds);
        end
    end

% save the data in simulation structure
    simulation.completeLP.runtime(exp) = time_completeLP;
    simulation.completeLP.cost(exp) = cost_completeLP;
    simulation.completeLP.cellCapacity(exp) = unrounded_s_completeLP;
    simulation.completeLP.sigma(exp) = {sigma_completeLP};

% save the data in simulation structure
    simulation.interMILP.runtime(exp) = 0;
    simulation.interMILP.cost(exp) = 0;
    simulation.interMILP.no_inter_markings(exp) = 0;
    simulation.interMILP.sigma(exp) = cell(1,1);

% save the data in simulation structure
    simulation.interLP.runtime(exp) = 0;
    simulation.interLP.cost(exp) = 0;
    simulation.interLP.no_inter_markings(exp) = 0;
    simulation.interLP.sigma(exp) = {sigma_completeLP};

% if s* is integer then sigma* should be equal with sigmaM - go to the next experiment
% if s* is greater than 1 then s_bar = ceil(s*) and solve the problem with intermediary markings

    if norm(sol_completeLP.s-1)>10*tol
%% (3) MILP with intermediary markings (for collision avoidance) - m and sigma as variables
% s_bar unknown markings (integer variables) + final marking is known (m_s_bar = V'*x)
% s_bar unknown firing vectors (integer variables)
% cell capacity is equal with 1 in order to ensure collision avoidance
        s_bar = ceil(sol_completeLP.s);
        
        m = optimvar('m',np,s_bar+1,'LowerBound',0,'UpperBound',N_r,'Type','integer');
        sigma = optimvar('sigma',nt,s_bar,'LowerBound',0,'UpperBound',inf,'Type','integer');

        prob = optimproblem;
        prob.Objective = sum(sum(sigma + sigma*diag(s_bar-(1:s_bar)')));

        prob.Constraints.cons1 = m(:,2:s_bar+1) == m(:,1:s_bar) + C*sigma;
        prob.Constraints.cons2 = Post * sigma + m(:,1:s_bar) <= ones(np,s_bar);
        prob.Constraints.cons3 = m(:,end) == V'*x;
        prob.Constraints.cons4 = m(:,1) == m0;

        opt = optimoptions('intlinprog','Display','none');

        tic;
        [sol_interMILP,fval_interMILP,exitflag_interMILP,output_interMILP] = solve(prob,'Solver','intlinprog','Options',opt);
        time_interMILP = toc;

        flag_notFeasible_interMILP = 0;
        count_notFeasible_interMILP = 0;
        if exitflag_interMILP~=1 && exitflag_interMILP~=2
            fprintf('MILP problem with intermediary markings is not feasbile.\n');
            count_notFeasible_interMILP = count_notFeasible_interMILP + 1;
            flag_notFeasible_interMILP = 1;
%             return;
        else
            if norm(sol_interMILP.sigma-round(sol_interMILP.sigma)) < 10*tol
                sigma_interMILP = round(sol_interMILP.sigma);
                m_interMILP = round(sol_interMILP.m);

                cost_interMILP = sum(sum(sigma_interMILP));

% compute trajectories                
                rob_traj_interMILP = compute_trajectories(T,Pre,Post,m_interMILP,sigma_interMILP,x0);
% plot trajectories 
                plot_environment_obstacles(regions,obstacles,env_bounds,reg_edges,T.Vert); %partition
                title('MILP problem with intermediary markings');
                plot_trajectories(rob_traj_interMILP,rob_color);
% plot animation                
                plot_animation(T,Pre,Post,m_interMILP,sigma_interMILP,env_bounds);

            else
                fprintf('MILP problem with intermediary markings - sigma is not integer!\n');
            end
        end

% save the data in simulation structure
        if flag_notFeasible_interMILP
            simulation.interMILP.runtime(exp) = time_interMILP;
            simulation.interMILP.cost(exp) = Inf;
            simulation.interMILP.no_inter_markings(exp) = Inf;
            simulation.interMILP.sigma(exp) = cell(1,1);
        else
            simulation.interMILP.runtime(exp) = time_interMILP;
            simulation.interMILP.cost(exp) = cost_interMILP;
            simulation.interMILP.no_inter_markings(exp) = s_bar;
            simulation.interMILP.sigma(exp) = {sigma_interMILP};
        end

%% (4) LP with intermediate markings (for collision avoidance) - m and sigma as variables
% s_bar unknown markings (continuous variables) + final marking is known (m_s_bar = V'*x)
% s_bar unknown firing vectors (continuous variables)
% cell capacity is equal with 1 in order to ensure collision avoidance
        sigma.Type = 'continuous';
        m.Type = 'continuous';
        
        opt = optimoptions('intlinprog','Display','none');
        
        tic;
        [sol_interLP,fval_interLP,exitflag_interLP,output_interLP] = solve(prob,'Solver','intlinprog','Options',opt);
        time_interLP = toc;
        
        flag_notFeasible_interLP = 0;
        count_notFeasible_interLP = 0;
        if exitflag_interLP~=1 && exitflag_interLP~=2
            fprintf('relaxed LP problem with intermediary markings is not feasbile.\n');
            count_notFeasible_interLP = count_notFeasible_interLP + 1;
            flag_notFeasible_interLP = 1;
%             return;
        else
            if norm(sol_interLP.sigma-round(sol_interLP.sigma)) < 10*tol
                sigma_interLP = round(sol_interLP.sigma);
                m_interLP = round(sol_interLP.m);
        
                cost_interLP = sum(sum(sigma_interLP));

% compute trajectories                
                rob_traj_interLP = compute_trajectories(T,Pre,Post,m_interLP,sigma_interLP,x0);
% plot trajectories 
                plot_environment_obstacles(regions,obstacles,env_bounds,reg_edges,T.Vert); %partition
                title('relaxed LP problem with intermediary markings');
                plot_trajectories(rob_traj_interLP,rob_color);
% plot animation                
                plot_animation(T,Pre,Post,m_interLP,sigma_interLP,env_bounds);

            else
                fprintf('relaxed LP problem with intermediary markings - sigma is not integer!\n');
            end
        end
% save the data in simulation structure
        if flag_notFeasible_interLP
            simulation.interLP.runtime(exp) = time_interLP;
            simulation.interLP.cost(exp) = Inf;
            simulation.interLP.no_inter_markings(exp) = Inf;
            simulation.interLP.sigma(exp) = cell(1,1);
        else
            simulation.interLP.runtime(exp) = time_interLP;
            simulation.interLP.cost(exp) = cost_interLP;
            simulation.interLP.no_inter_markings(exp) = s_bar;
            simulation.interLP.sigma(exp) = {sigma_interLP};
        end
    else
        fprintf('Cell capacity is equal with 1 >>> no intermediary markings are needed.\n');
    end
end
