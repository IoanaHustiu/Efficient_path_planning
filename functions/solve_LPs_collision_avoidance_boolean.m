function [simulation,flag] = solve_LPs_collision_avoidance_boolean(Post,Pre,At,bt,m0,T,flag_ILP,plot_animation)
flag = 1;
nplaces = size(Post,1);
ntrans = size(Post,2);
nR = sum(m0);
env_limit = size(T.map2D);
np = size(At,2);
V = zeros(np,nplaces);
idx = sub2ind([np,nplaces],(1:np)',T.props(:));
V(idx) = 1;

if flag_ILP
    [MILP1,ILP1] = solve_mILP_boolean(Post,Pre,V,At,bt,m0,flag_ILP);
    if ILP1.exitflag~=1
        if MILP1.exitflag==10 && ILP1.exitflag==10
            simulation.ILP.timelimit = 0;
            simulation.MILP1.timelimit = 0;
            simulation.ILP.exitflag = ILP1.exitflag;
            flag = 0;
            return;
        else
            fprintf('ILP >>> exitflag ILP = %i.\n',ILP1.exitflag);
            simulation.ILP.timelimit = 0;
            simulation.completeILP.runtime = ILP1.runtime;
            simulation.completeILP.cost = 0;
            simulation.completeILP.sol = 0;
            simulation.completeILP.cellCapacity = 0;
            simulation.ILP.exitflag = ILP1.exitflag;
        end
    else
        simulation.ILP.timelimit = 1;
        simulation.completeILP.runtime = ILP1.runtime;
        simulation.completeILP.cost = sum(ILP1.sol(nplaces+1:nplaces+ntrans));
        simulation.completeILP.sol = ILP1.sol;
        simulation.completeILP.cellCapacity = ILP1.sol(end);
        simulation.ILP.exitflag = ILP1.exitflag;
        s_ILP = ceil(ILP1.sol(end)); % cell capacity
    end
else
    MILP1 = solve_mILP_boolean(Post,Pre,V,At,bt,m0,flag_ILP);
end

if MILP1.exitflag~=1
    if MILP1.exitflag==10
        simulation.ILP.timelimit = 0;
        simulation.MILP1.timelimit = 0;
        simulation.completeLP.exitflag = MILP1.exitflag;
        flag = 0;
        return;
    else
        fprintf('LP1 >>> exitflag LP1 = %i.\n',MILP1.exitflag);
        simulation.MILP1.timelimit = 0;
        simulation.completeLP.exitflag = MILP1.exitflag;
        flag = 0;
        return;
    end
else
    x = MILP1.sol(nplaces+ntrans+1:nplaces+ntrans+size(At,2));
    s = MILP1.sol(end); % cell capacity
end

% save the data in simulation structure
simulation.completeLP.runtime = MILP1.runtime;
simulation.completeLP.cost = sum(MILP1.sol(nplaces+1:nplaces+ntrans));
simulation.completeLP.sol = MILP1.sol;
simulation.completeLP.cellCapacity = MILP1.sol(end);
simulation.completeLP.exitflag = MILP1.exitflag;

%check if the cell capacity is equal with 1
if norm(s - 1) < 1000*eps && norm(x - round(x)) < 1000*eps %cell capacity is 1 and collision avoidance is not required
    num_intermediate = 1;
    fprintf('Cell capacity is equal with 1 and x is integer. Second problem will not be solved!\n');

    fprintf("\nCollision avoidance is already imposed. The second problem will not be solved.\n");
    % save the data in simulation structure
    simulation.interLP.runtime = 0;
    simulation.interLP.cost = 0;
    simulation.interLP.sol = 0;
    simulation.interLP.no_inter_markings = 0;

    if flag_ILP
        % save the data in simulation structure
        simulation.ILPIM.runtime = 0;
        simulation.ILPIM.cost = 0; simulation.completeLP.cost
    end

    fprintf(1,'\n============================================================');
    fprintf(1,'\n    Solution for LP1 obtained using intlinprog solver ');
    fprintf(1,'\n============================================================');
    fprintf(1,'\nSolution for LP1 found in %f [secs].',MILP1.runtime);
    fprintf(1,'\nOptimal solution LP1 = %s',num2str(MILP1.cost));
    fprintf(1,'\nInfinite norm for LP1 of Post * sigma = %s',num2str(MILP1.sol(end)));
else
    if plot_animation
        fprintf(1,'\n============================================================');
        fprintf(1,'\n    Solution for LP1 obtained using intlinprog solver ');
        fprintf(1,'\n============================================================');
        fprintf(1,'\nSolution for LP1 found in %f [secs].',MILP1.runtime);
        fprintf(1,'\nOptimal solution LP1 = %s',num2str(MILP1.cost));
        fprintf(1,'\nInfinite norm for LP1 of Post * sigma = %s',num2str(MILP1.sol(end)));
        %solve the problem with intermediate markings
        fprintf(1,'\n\n=============================================');
        fprintf(1,'\n           Start solving second problem...');
        fprintf(1,'\n=============================================');
    end

    clear Aeq1; clear beq1;
    clear Aineq1; clear bineq1;

    num_intermediate = ceil(s);

    flag_num = 1;

    while flag_num
        fprintf(1,'\nSolve second problem with %d intermediate markings\n',num_intermediate);

        if flag_ILP && simulation.ILP.timelimit
            [LPIM,ILPIM] = solve_LP_IM_boolean(num_intermediate,Post-Pre,Post,At,bt,V,x,m0,flag_ILP);

            if ILPIM.exitflag~=1
                if ILPIM.exitflag == 0
                    if num_intermediate < nR
                        num_intermediate = num_intermediate + 1;
                    end
                elseif ILPIM.exitflag == 10 && LPIM.exitflag == 10
                    flag = 0;
                    flag_num = 0;
                    return;
                elseif ILPIM.exitflag == 3600
                    simulation.ILPIM.timelimit = 0;
                end
            else
                if ILP1.runtime + ILPIM.runtime < 3600
                    flag_num = 0;
                    simulation.ILPIM.success = 1;
                    % save the data in simulation structure
                    simulation.ILPIM.runtime = ILPIM.runtime;
                    simulation.ILPIM.cost = 0;
                    simulation.ILPIM.sol = ILPIM.sol;
                else
                    simulation.ILPIM.success = 0;
                    % return;
                end
            end
        else
            LPIM = solve_LP_IM_boolean(num_intermediate,Post-Pre,Post,At,bt,V,x,m0,flag_ILP);

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
                        flag = 0;
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
    end

    % save the data in simulation structure
    simulation.interLP.runtime = LPIM.runtime;
    simulation.interLP.cost = 0;
    simulation.interLP.sol = LPIM.sol;
    simulation.interLP.no_inter_markings = num_intermediate;

    for i = 1 : num_intermediate
        simulation.interLP.cost = simulation.interLP.cost + sum(LPIM.sol((i-1)*(nplaces+ntrans) + nplaces + 1 : i*(nplaces+ntrans)));
        if flag_ILP
            simulation.ILPIM.cost = simulation.ILPIM.cost + sum(ILPIM.sol((i-1)*(nplaces+ntrans) + nplaces + 1 : i*(nplaces+ntrans)));
        end
    end

    if flag_ILP
        if simulation.interLP.cost ~= simulation.ILPIM.cost
            fprintf(1,'===============\n');
            fprintf(1,'(8) - MILP (m_i, sigma_i, x are unknowns) cost: %i\n',simulation.interLP.cost);
            fprintf(1,'(8) - ILP (m_i, sigma_i and x are unknowns) cost: %i\n',simulation.ILPIM.cost);
            fprintf(1,'===============\n');
        else
            fprintf(1,'===============\n');
            fprintf(1,'Optimal cost is obtained! (8) - ILP cost: %i, (8) - MILP cost %i\n',simulation.ILPIM.cost,simulation.interLP.cost);
            fprintf(1,'===============\n');
        end
    end
end

fprintf(1,'\nSolution MILP found in %f [secs].',simulation.interLP.runtime);
fprintf(1,'\nOptimal value MILP = %s\n',num2str(simulation.interLP.cost));
fprintf(1,'\nTotal time for solution LP1+MILP: %f [sec].',simulation.completeLP.runtime+simulation.interLP.runtime);

if flag_ILP
    fprintf(1,'\nTotal time for solution ILP+ILP with intermediary markings: %f [sec].',simulation.completeILP.runtime+simulation.ILPIM.runtime);
end

if (plot_animation)
    fprintf(1,'\n============================================================');
    fprintf(1,'\n    Solution for MILP obtained using intlinprog solver ');
    fprintf(1,'\n============================================================\n');
    % fprintf(1,'\nSolution for the MILP found in %f [secs].',MILP1.runtime);
    % fprintf(1,'\nOptimal solution for the MILP = %s',num2str(MILP1.cost));
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
        % figure;
        plot_environment(T.rem_cell(find(current_marking)),T.rem_cell(T.props),T.map2D,env_limit,T.Vert); 
        title(sprintf('Efficient path planning (LP1 + MILP): iteration %d',i));
        hold on;

        if num_intermediate == 1
            current_sigma = round(MILP1.sol(nplaces+1:nplaces+ntrans));
        else
            current_sigma = round(LPIM.sol((nplaces+ntrans)*(i-1)+nplaces+1:(nplaces+ntrans)*(i-1)+nplaces+ntrans));
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

        for j = 1 : length(Rob_places)
            traj = Rob_places{j};
            for k = 1 : length(traj)-1
                plot([T.centr{T.rem_cell(traj(k))}(1) T.centr{T.rem_cell(traj(k+1))}(1)],...
                    [T.centr{T.rem_cell(traj(k))}(2) T.centr{T.rem_cell(traj(k+1))}(2)],'-','LineWidth',1,'Color',rob_color(j,:));
            end
            if (length(traj) > 1)
               init = T.centr{T.rem_cell(traj(1))};
               fill([init(1)-0.5 init(1)+0.5 init(1)+0.5 init(1)-0.5],[init(2)-0.5 init(2)-0.5 init(2)+0.5 init(2)+0.5],rob_color(j,:),'EdgeColor','none');
               init = T.centr{T.rem_cell(traj(end))};
               fill([init(1)-0.5 init(1)+0.5 init(1)+0.5 init(1)-0.5],[init(2)-0.5 init(2)-0.5 init(2)+0.5 init(2)+0.5],rob_color(j,:),'EdgeColor','none');
            end
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
        end
        saveas(gcf,sprintf('smart_plant_eff_it%d.fig', i));
    end

    if flag_ILP && simulation.ILP.timelimit
        rob_color = hsv(sum(m0));

        for i = 1 : num_intermediate
            if (i == 1)
                current_marking = m0;
            else
                current_marking = round(ILPIM.sol((nplaces+ntrans)*(i-2)+1:(nplaces+ntrans)*(i-2)+nplaces));
            end
            %figure;
            plot_environment(T.rem_cell(find(current_marking)),T.rem_cell(T.props),T.map2D,env_limit,T.Vert);
            title(sprintf('ILP formulation: iteration %d',i));
            hold on;

            if num_intermediate == 1
                current_sigma = round(ILP1.sol(nplaces+1:nplaces+ntrans));
            else
                current_sigma = round(ILPIM.sol((nplaces+ntrans)*(i-1)+nplaces+1:(nplaces+ntrans)*(i-1)+nplaces+ntrans));
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

            for j = 1 : length(Rob_places)
                traj = Rob_places{j};
                for k = 1 : length(traj)-1
                    plot([T.centr{T.rem_cell(traj(k))}(1) T.centr{T.rem_cell(traj(k+1))}(1)],...
                        [T.centr{T.rem_cell(traj(k))}(2) T.centr{T.rem_cell(traj(k+1))}(2)],'-','LineWidth',1,'Color',rob_color(j,:));
                end
                if (length(traj) > 1)
                    init = T.centr{T.rem_cell(traj(1))};
                    fill([init(1)-0.5 init(1)+0.5 init(1)+0.5 init(1)-0.5],[init(2)-0.5 init(2)-0.5 init(2)+0.5 init(2)+0.5],rob_color(j,:),'EdgeColor','none');
                    init = T.centr{T.rem_cell(traj(end))};
                    fill([init(1)-0.5 init(1)+0.5 init(1)+0.5 init(1)-0.5],[init(2)-0.5 init(2)-0.5 init(2)+0.5 init(2)+0.5],rob_color(j,:),'EdgeColor','none');
                end
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
            end
            saveas(gcf,sprintf('smart_plant_ILP_it%d.fig', i));
        end
    end
end