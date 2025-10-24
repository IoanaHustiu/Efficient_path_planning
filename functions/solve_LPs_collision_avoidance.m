function [simulation,flag] = solve_LPs_collision_avoidance(Post,Pre,mf,m0,T,flag_ILP,plot_animation)
flag = 1;
nplaces = size(Post,1);
ntrans = size(Post,2);
nR = sum(m0);
np = nR;
env_limit = size(T.map2D);

fprintf(1,'\nSolve first problem\n');
if flag_ILP
    [MILP1,ILP1] = solve_mILPr(Post,Pre,mf,m0,flag_ILP);

    if ILP1.exitflag~=1
        if MILP1.exitflag==10 && ILP1.exitflag==10
            flag = 0;
            return;
        elseif ILP1.exitflag~=3600
            fprintf('Unfeasible problem - exitflag ILP ~= 1.\n');
            flag = 0;
            return;
        else
            flagILP = 0; %time limit has been exceeded
            simulation.ILP.timelimit = 0;
            return;
        end
    else
        simulation.completeILP.runtime = ILP1.runtime;
        simulation.completeILP.cost = sum(ILP1.sol(1:ntrans));
        simulation.completeILP.sol = ILP1.sol(1:ntrans);
        simulation.completeILP.cellCapacity = ILP1.sol(end);
        s_ILP = ceil(ILP1.sol(end)); % cell capacity
    end

else
    MILP1 = solve_mILPr(Post,Pre,mf,m0,flag_ILP);
end

if MILP1.exitflag~=1
    if MILP1.exitflag==10
        flag = 0;
        return;
    elseif MILP1.exitflag~=3600
        fprintf('Unfeasible problem - exitflag LP1 ~= 1.\n');
        flag = 0;
        return;
    else
        flag = 0; %time limit has been exceeded
        simulation.MILP1.timelimit = 0;
        return;
    end
else
    s = ceil(MILP1.sol(end)); % cell capacity
end

% save the data in simulation structure

simulation.completeLP.runtime = MILP1.runtime;
simulation.completeLP.cost = sum(MILP1.sol(1:ntrans));
simulation.completeLP.sol = MILP1.sol(1:ntrans);
simulation.completeLP.cellCapacity = MILP1.sol(end);

%check if the cell capacity is equal with 1
if norm(s - 1) < 1000*eps %cell capacity is 1 and collision avoidance is not required
    fprintf('Cell capacity is equal with 1. No collision avoidance is necessary!\n');

    fprintf(1,'\n============================================================');
    fprintf(1,'\n    Solution for LP1/IP1 obtained using intlinprog solver ');
        fprintf(1,'\n--------------------------------------------------------');
    fprintf(1,'\nSolution for LP1 found in %f [secs].',MILP1.runtime);
    fprintf(1,'\nOptimal solution for LP1 = %s',num2str(MILP1.cost));
    fprintf(1,'\nInfinite norm for LP1 of Post * sigma = %s',num2str(MILP1.sol(end)));
    if flag_ILP
        fprintf(1,'\n--------------------------------------------------------');
        fprintf(1,'\nSolution for IP1 found in %f [secs].',ILP1.runtime);
        fprintf(1,'\nOptimal solution for IP1 = %s',num2str(ILP1.cost));
        fprintf(1,'\nInfinite norm for IP1 of Post * sigma = %s',num2str(ILP1.sol(end)));
    end
end

if s > 1
    fprintf(1,'\n============================================================');
    fprintf(1,'\n    Solution for LP1 obtained using intlinprog solver ');
    fprintf(1,'\n============================================================');
    fprintf(1,'\nSolution for LP1 found in %f [secs].',MILP1.runtime);
    fprintf(1,'\nOptimal solution for LP1 = %s',num2str(MILP1.cost));
    fprintf(1,'\nInfinite norm for LP1 of Post * sigma = %s',num2str(MILP1.sol(end)));
    if flag_ILP
        fprintf(1,'\n--------------------------------------------------------');
        fprintf(1,'\nSolution for IP1 found in %f [secs].',ILP1.runtime);
        fprintf(1,'\nOptimal solution for IP1 = %s',num2str(ILP1.cost));
        fprintf(1,'\nInfinite norm for IP1 of Post * sigma = %s',num2str(ILP1.sol(end)));
    end
    %solve the problem with intermediate markings
    fprintf(1,'\n\n=============================================');
    fprintf(1,'\n           Start solving second problem...');
    fprintf(1,'\n=============================================');
end

clear Aeq1; clear beq1;
clear Aineq1; clear bineq1;

num_intermediate = s;
if num_intermediate > 1
    flag_num = 1;

    while flag_num
        fprintf(1,'\nSolve second problem with %d intermediate markings\n',num_intermediate);

        if flag_ILP
            [LPIM,ILPIM] = solve_LP_IMr(num_intermediate,Post-Pre,Post,m0,mf,flag_ILP);

            if ILPIM.exitflag~=1
                if ILPIM.exitflag == 10 && LPIM.exitflag == 10
                    flag = 0;
                    flag_num = 0;
                    return;
                elseif ILPIM.exitflag == 3600
                    simulation.ILPIM.timelimit = 0;
                end
            else
                if ILP1.runtime + ILPIM.runtime < 3600
                    simulation.ILPIM.success = 1;
                    % save the data in simulation structure
                    simulation.ILPIM.runtime = ILPIM.runtime;
                    simulation.ILPIM.cost = 0;
                    simulation.ILPIM.sol = ILPIM.sol;
                else
                    simulation.ILPIM.success = 0;
                    return;
                end
            end
        else
            LPIM = solve_LP_IMr(num_intermediate,Post-Pre,Post,m0,mf,flag_ILP);
        end

        if LPIM.exitflag~=1
            if LPIM.exitflag == 10
                flag = 0;
                flag_num = 0;
                return;
            elseif LPIM.exitflag == 3600
                flag = 0; %time limit has been exceeded
                simulation.timelimit = 0;
            else
                num_intermediate = num_intermediate + 1;
                if num_intermediate > nR
                    fprintf('Number of intermediary markings is greater than the number of robots.\n');
                    fprintf('Solution not found!\n');
                    flag_num = 0;
                    return;
                end
            end
        else
            flag_num = 0;
            if MILP1.runtime + LPIM.runtime < 3600
                simulation.success = 1;
            else
                simulation.success = 0;
                return;
            end
        end
    end

    % save the data in simulation structure
    simulation.interLP.runtime = LPIM.runtime;
    simulation.interLP.cost = 0;
    simulation.interLP.sol = LPIM.sol;
    simulation.interLP.no_inter_markings = num_intermediate;

    for i = 1 : num_intermediate
        simulation.interLP.cost = simulation.interLP.cost + sum(LPIM.sol((i-1)*(nplaces+ntrans) + 1 + nplaces:i*(nplaces+ntrans)));
        if flag_ILP
            simulation.ILPIM.cost = simulation.ILPIM.cost + sum(ILPIM.sol((i-1)*(nplaces+ntrans) + 1 + nplaces:i*(nplaces+ntrans)));
        end
    end
else
    fprintf("\nCollision avoidance is already imposed. The second problem will not be solved.\n");
    % save the data in simulation structure
    simulation.interLP.runtime = 0;
    simulation.interLP.cost = 0;
    simulation.interLP.sol = 0;
    simulation.interLP.no_inter_markings = 0;

    if flag_ILP
        % save the data in simulation structure
        simulation.ILPIM.runtime = 0;
        simulation.ILPIM.cost = 0;
    end
end

if (num_intermediate > 1)
    fprintf(1,'\nSolution LP2 found in %f [secs].',simulation.interLP.runtime);
    fprintf(1,'\nOptimal value LP2 = %s\n',num2str(simulation.interLP.cost));
    fprintf(1,'Total time for solution LP1+LP2: %f [sec].\n',simulation.completeLP.runtime+simulation.interLP.runtime);
    if flag_ILP
        fprintf(1,'\nSolution IP2 found in %f [secs].',simulation.ILPIM.runtime);
        fprintf(1,'\nOptimal value IP2 = %s\n',num2str(simulation.ILPIM.cost));
        fprintf(1,'Total time for solution ILP1+ILP2: %f [sec].\n',simulation.completeILP.runtime+simulation.ILPIM.runtime);
    end
else
    fprintf(1,'\nTotal time for solution LP1: %f [sec].',simulation.completeLP.runtime+simulation.interLP.runtime);
    if flag_ILP
        fprintf(1,'\nTotal time for solution ILP1: %f [sec].',simulation.completeILP.runtime+simulation.ILPIM.runtime);
    end
end

if plot_animation
    %rob_color={'b','g','c','m','k','y'};
    rob_color = hsv(sum(m0));

    for i = 1 : num_intermediate
        if (i == 1)
            current_marking = m0;
        else
            current_marking = round(LPIM.sol((nplaces+ntrans)*(i-2)+1:(nplaces+ntrans)*(i-2)+nplaces));
        end
        %figure;
        %        plot_environment(T.rem_cells(find(current_marking)),T.rem_cells(T.props),T.map2D,env_limit,T.Vert);
        % --- Extract indices of start and goal cells
        start_idx = T.rem_cells(find(current_marking));  % starts: places with tokens
        goal_idx  = T.rem_cells(T.props);               % goals: places in T.props

        % --- Ensure same number of start and goal points
        n_r = min(numel(start_idx), numel(goal_idx));
        if n_r == 0
            warning('No start/goal pairs to plot.');
            return;
        end
        start_idx = start_idx(1:n_r);
        goal_idx  = goal_idx(1:n_r);

        % --- Convert linear indices to (row,col) then to 0-based Cartesian [x y]
        [H, W] = size(T.map2D);
        [start_r, start_c] = ind2sub([H, W], start_idx(1:n_r));
        [goal_r,  goal_c ] = ind2sub([H, W], goal_idx (1:n_r));

        % correct centers (x = c - 0.5  →  +0.5 for Cartesian surface alignment)
        selectedStart = [start_c - 0.5 , start_r - 0.5 ];
        selectedFin   = [goal_c  - 0.5 , goal_r  - 0.5 ];

        % --- Plot environment with starts & goals
        %plot_environment_new(selectedPts, T.map2D, T);
        plot_environment_new_SG(selectedStart, selectedFin, T.map2D, T);

        title(sprintf('Efficient path planning (LP1 + LP2): iteration %d',i));
        hold on;

        if num_intermediate == 1
            current_sigma = round(MILP1.sol(1:ntrans));
        else
            base = (nplaces+ntrans)*(i-1);
            current_sigma = round(LPIM.sol(base + nplaces + (1:ntrans)));
        end

        try
            [feasible_sigma, Rob_places, ~, ~] = sigma2trajectories(Pre,Post,current_marking,current_sigma,find(current_marking));
            if ~feasible_sigma
                error('something wrong happened!');
            end
        catch
            warning('something wrong happened!');
            flag = 0;
            return;
        end
        nR = numel(Rob_places);
        if nR == 0
            continue;
        end
        % Asegura suficientes colores
        if size(rob_color,1) < nR
            rob_color = hsv(nR);
        end
        % Dibuja cada trayectoria
        for j = 1:nR
            traj = Rob_places{j};
            XY = places2xy(traj, T);        % Nx2 con [x y] de cada lugar

            % Línea de la trayectoria
            plot(XY(:,1), XY(:,2), '-', 'LineWidth', 1.5, 'Color', rob_color(j,:));

            % Marcadores inicio/fin
            c_start = XY(1,:);
            c_goal  = XY(end,:);
            scatter(c_start(1), c_start(2), 36, rob_color(j,:), '^', 'filled', 'MarkerEdgeColor',rob_color(j,:));  % inicio
            scatter(c_goal(1),  c_goal(2),  36, rob_color(j,:), 'd', 'filled', 'MarkerEdgeColor',rob_color(j,:));  % fin
        end
        %for j = 1 : length(Rob_places)
        %    traj = Rob_places{j};
        %    for k = 1 : length(traj)-1
        %        plot([T.centr{T.rem_cells(traj(k))}(1) T.centr{T.rem_cells(traj(k+1))}(1)],...
        %            [T.centr{T.rem_cells(traj(k))}(2) T.centr{T.rem_cells(traj(k+1))}(2)],'-','LineWidth',1,'Color',rob_color(j,:));
        %    end
        %    if (length(traj) > 1)
        %       init = T.centr{T.rem_cells(traj(1))};
        %       fill([init(1)-0.5 init(1)+0.5 init(1)+0.5 init(1)-0.5],[init(2)-0.5 init(2)-0.5 init(2)+0.5 init(2)+0.5],rob_color(j,:),'EdgeColor','none');
        %       init = T.centr{T.rem_cells(traj(end))};
        %       fill([init(1)-0.5 init(1)+0.5 init(1)+0.5 init(1)-0.5],[init(2)-0.5 init(2)-0.5 init(2)+0.5 init(2)+0.5],rob_color(j,:),'EdgeColor','none');
        %    end
            %final destinations already reached
            % intersectt = intersect(traj(end),find(mf));
            % if ~isempty(intersectt)
            %     reached = T.centr{T.rem_cell(intersectt)};
            %     %text(reached(1), reached(2), 'X', 'HorizontalAlignment','center','FontSize',12);
            %     c = reached;
            %     theta = linspace(0,2*pi,64);
            %     r = 1;
            %     x = c(1) + r*cos(theta);
            %     y = c(2) + r*sin(theta);
            %     plot(x, y, 'g-', 'LineWidth', 2);
            % end
        %end
    end

    % if flag_ILP
    %     rob_color = hsv(sum(m0));
    % 
    %     for i = 1 : num_intermediate
    %         if (i == 1)
    %             current_marking = m0;
    %         else
    %             current_marking = round(ILPIM.sol((nplaces+ntrans)*(i-2)+1:(nplaces+ntrans)*(i-2)+nplaces));
    %         end
    %         %figure;
    %         plot_environment(T.rem_cells(find(current_marking)),T.rem_cells(T.props),T.map2D,env_limit,T.Vert);
    %         title(sprintf('ILP formulation: iteration %d',i));
    %         hold on;
    % 
    %         if num_intermediate == 1
    %             current_sigma = round(ILP1.sol(1:ntrans));
    %         else
    %             current_sigma = round(ILPIM.sol((nplaces+ntrans)*(i-1) + nplaces + 1 : (nplaces+ntrans)*(i-1) + nplaces + ntrans));
    %         end
    % 
    %         try
    %             [feasible_sigma, Rob_places, ~, ~] = sigma2trajectories(Pre,Post,current_marking,current_sigma,find(current_marking));
    %             if ~feasible_sigma
    %                 error('something wrong happened!');
    %             end
    %         catch
    %             warning('something wrong happened!');
    %             flag = 0;
    %             return;
    %         end
    % 
    %         for j = 1 : length(Rob_places)
    %             traj = Rob_places{j};
    %             for k = 1 : length(traj)-1
    %                 plot([T.centr{T.rem_cells(traj(k))}(1) T.centr{T.rem_cells(traj(k+1))}(1)],...
    %                     [T.centr{T.rem_cells(traj(k))}(2) T.centr{T.rem_cells(traj(k+1))}(2)],'-','LineWidth',1,'Color',rob_color(j,:));
    %             end
    %             if (length(traj) > 1)
    %                 init = T.centr{T.rem_cells(traj(1))};
    %                 fill([init(1)-0.5 init(1)+0.5 init(1)+0.5 init(1)-0.5],[init(2)-0.5 init(2)-0.5 init(2)+0.5 init(2)+0.5],rob_color(j,:),'EdgeColor','none');
    %                 init = T.centr{T.rem_cells(traj(end))};
    %                 fill([init(1)-0.5 init(1)+0.5 init(1)+0.5 init(1)-0.5],[init(2)-0.5 init(2)-0.5 init(2)+0.5 init(2)+0.5],rob_color(j,:),'EdgeColor','none');
    %             end
    %             %final destinations already reached
    %             % intersectt = intersect(traj(end),find(mf));
    %             % if ~isempty(intersectt)
    %             %     reached = T.centr{T.rem_cell(intersectt)};
    %             %     %text(reached(1), reached(2), 'X', 'HorizontalAlignment','center','FontSize',12);
    %             %     c = reached;
    %             %     theta = linspace(0,2*pi,64);
    %             %     r = 1;
    %             %     x = c(1) + r*cos(theta);
    %             %     y = c(2) + r*sin(theta);
    %             %     plot(x, y, 'g-', 'LineWidth', 2);
    %             % end
    %         end
    %     end
    % end

end