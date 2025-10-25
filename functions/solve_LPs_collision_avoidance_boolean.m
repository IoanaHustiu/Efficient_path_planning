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
    [LP1,ILP1] = solve_mILP_boolean(Post,Pre,V,At,bt,m0,flag_ILP);
    if ILP1.exitflag~=1
        if LP1.exitflag==10 && ILP1.exitflag==10
            simulation.ILP1.timelimit = 0;
            simulation.ILP1.runtime = 0;
            simulation.ILP1.cost = 0;
            simulation.ILP1.sol = 0;
            simulation.ILP1.cellCapacity = 0;
            simulation.ILP1.exitflag = ILP1.exitflag;

            simulation.LP1.timelimit = 0;
            simulation.LP1.runtime = 0;
            simulation.LP1.cost = 0;
            simulation.LP1.sol = 0;
            simulation.LP1.cellCapacity = 0;
            simulation.LP1.exitflag = LP1.exitflag;
            flag = 0;
            return;
        elseif ILP1.exitflag == 3600
            simulation.ILP1.timelimit = 0;
            simulation.ILP1.runtime = ILP1.runtime;
            simulation.ILP1.exitflag = ILP1.exitflag;
        else
            fprintf('ILP >>> exitflag ILP = %i.\n',ILP1.exitflag);
            simulation.ILP1.timelimit = 0;
            simulation.ILP1.runtime = ILP1.runtime;
            simulation.ILP1.cost = 0;
            simulation.ILP1.sol = 0;
            simulation.ILP1.cellCapacity = 0;
            simulation.ILP1.exitflag = ILP1.exitflag;
        end
    else
        simulation.ILP1.timelimit = 1;
        simulation.ILP1.runtime = ILP1.runtime;
        simulation.ILP1.cost = sum(ILP1.sol(nplaces+1:nplaces+ntrans));
        simulation.ILP1.sol = ILP1.sol;
        simulation.ILP1.cellCapacity = ILP1.sol(end);
        simulation.ILP1.exitflag = ILP1.exitflag;
    end
else
    LP1 = solve_mILP_boolean(Post,Pre,V,At,bt,m0,flag_ILP);
end

if LP1.exitflag~=1
    if LP1.exitflag==10
        simulation.LP1.timelimit = 0;
        simulation.LP1.runtime = 0;
        simulation.LP1.exitflag = LP1.exitflag;
        simulation.LP1.cost = 0;
        simulation.LP1.sol = 0;
        simulation.LP1.cellCapacity = 0;
    elseif LP1.exitflag == 3600
        simulation.LP1.timelimit = 0;
        simulation.LP1.runtime = LP1.runtime;
        simulation.LP1.exitflag = LP1.exitflag;
    else
        fprintf('LP1 >>> exitflag LP1 = %i.\n',LP1.exitflag);
        simulation.LP1.timelimit = 0;
        simulation.LP1.runtime = LP1.runtime;
        simulation.LP1.exitflag = LP1.exitflag;
        simulation.LP1.cost = 0;
        simulation.LP1.sol = 0;
        simulation.LP1.cellCapacity = 0;
    end
    flag = 0;
    return;
else
    x = LP1.sol(nplaces+ntrans+1:nplaces+ntrans+np);
    s = LP1.sol(end); % cell capacity

    % save the data in simulation structure
    simulation.LP1.timelimit = 1;
    simulation.LP1.runtime = LP1.runtime;
    simulation.LP1.cost = sum(LP1.sol(nplaces+1:nplaces+ntrans));
    simulation.LP1.sol = LP1.sol;
    simulation.LP1.cellCapacity = LP1.sol(end);
    simulation.LP1.exitflag = LP1.exitflag;
end

%check if the cell capacity is equal with 1
if norm(s - 1) < 1000*eps && norm(x - round(x)) < 1000*eps %cell capacity is 1 and x vector is integer
    num_intermediate = 1;
    fprintf('Cell capacity is equal with 1 and x is integer. Second problem will not be solved!\n');

    fprintf("\nCollision avoidance is already imposed. The second problem will not be solved.\n");
    % save the data in simulation structure
    simulation.MILP.runtime = 0;
    simulation.MILP.cost = 0;
    simulation.MILP.sol = 0;
    simulation.MILP.num_intermediate = 0;

    if flag_ILP
        % save the data in simulation structure
        simulation.ILP2.runtime = 0;
        simulation.ILP2.cost = 0;
    end

    fprintf(1,'\n============================================================');
    fprintf(1,'\n    Solution for LP1 obtained using intlinprog solver ');
    fprintf(1,'\n============================================================');
    fprintf(1,'\nSolution for LP1 found in %f [secs].',LP1.runtime);
    fprintf(1,'\nOptimal solution LP1 = %s',num2str(LP1.cost));
    fprintf(1,'\nInfinite norm for LP1 of Post * sigma = %s',num2str(LP1.sol(end)));
else
    if plot_animation
        fprintf(1,'\n============================================================');
        fprintf(1,'\n    Solution for LP1 obtained using intlinprog solver ');
        fprintf(1,'\n============================================================');
        fprintf(1,'\nSolution for LP1 found in %f [secs].',LP1.runtime);
        fprintf(1,'\nOptimal solution LP1 = %s',num2str(LP1.cost));
        fprintf(1,'\nInfinite norm for LP1 of Post * sigma = %s',num2str(LP1.sol(end)));
        %solve the problem with intermediate markings
        fprintf(1,'\n\n=============================================');
        fprintf(1,'\n           Start solving second problem...');
        fprintf(1,'\n=============================================');
    end

    clear Aeq1; clear beq1;
    clear Aineq1; clear bineq1;

    num_intermediate = ceil(s);

    flag_num = 1;
    simulation.ILP2.runtime_all = 0;
    simulation.MILP.runtime_all = 0;
    simulation.ILP2.cost = 0;
    simulation.MILP.cost = 0;

    while flag_num
        fprintf(1,'\nSolve second problem with %d intermediate markings\n',num_intermediate);

        if flag_ILP
            [MILP,ILP2] = solve_LP_IM_boolean(num_intermediate,Post-Pre,Post,At,bt,V,m0,flag_ILP);

            if ILP2.exitflag~=1
                if ILP2.exitflag == 10 && MILP.exitflag == 10
                    simulation.ILP2.success = 0;
                    simulation.MILP.success = 0;
                    simulation.ILP2.exitflag = ILP2.exitflag;
                    simulation.MILP.exitflag = MILP.exitflag;
                    simulation.ILP2.num_intermediate = num_intermediate;
                    simulation.MILP.num_intermediate = num_intermediate;

                    flag = 0;
                    return;
                elseif ILP2.exitflag == 3600
                    simulation.ILP2.timelimit = 3600;
                    simulation.ILP2.success = 0;
                    simulation.ILP2.runtime = ILP2.runtime;
                    simulation.ILP2.exitflag = ILP2.exitflag;
                    simulation.ILP2.num_intermediate = num_intermediate;
                else
                    if num_intermediate < nR
                        num_intermediate = num_intermediate + 1;
                    else
                        simulation.ILP2.num_intermediate = num_intermediate;
                        simulation.MILP.num_intermediate = num_intermediate;
                        simulation.ILP2.success = 0;
                        simulation.MILP.success = 0;
                        simulation.ILP2.exitflag = ILP2.exitflag;
                        simulation.MILP.exitflag = MILP.exitflag;

                        flag = 0;
                        return;
                    end
                    simulation.ILP2.runtime_all = simulation.ILP2.runtime_all + ILP2.runtime;
                    simulation.MILP.runtime_all =  simulation.MILP.runtime_all + MILP.runtime;
                end
            else
                if ILP1.runtime + ILP2.runtime < 3600
                    flag_num = 0;
                    simulation.ILP2.success = 1;
                    simulation.ILP2.runtime = ILP2.runtime;
                    simulation.ILP2.cost = 0;
                    simulation.ILP2.sol = ILP2.sol;
                    simulation.ILP2.exitflag = ILP2.exitflag;
                    simulation.ILP2.num_intermediate = num_intermediate;
                else
                    flag_num = 0;
                    simulation.ILP2.success = 3600;
                    simulation.ILP2.runtime = ILP2.runtime;
                    simulation.ILP2.exitflag = ILP2.exitflag;
                    simulation.ILP2.cost = 0;
                    simulation.ILP2.sol = 0;
                    simulation.ILP2.num_intermediate = num_intermediate;
                    % return;
                end
            end

            if MILP.exitflag~=1
                if MILP.exitflag == 10
                    simulation.MILP.success = 0;
                    simulation.MILP.exitflag = MILP.exitflag;
                    simulation.MILP.num_intermediate = num_intermediate;
                    flag = 0;
                    return;
                elseif MILP.exitflag == 3600
                    simulation.MILP.success = 0;
                    simulation.MILP.exitflag = MILP.exitflag;
                    simulation.MILP.runtime = MILP.runtime;
                    simulation.MILP.num_intermediate = num_intermediate;
                    flag = 0; %time limit has been exceeded
                    return;
                end
            else
                % save the data in simulation structure
                simulation.MILP.success = 1;
                simulation.MILP.runtime = MILP.runtime;
                simulation.MILP.cost = 0;
                simulation.MILP.sol = MILP.sol;
                simulation.MILP.num_intermediate = num_intermediate;
                simulation.MILP.exitflag = MILP.exitflag;
                flag_num = 0;
            end
        else
            MILP = solve_LP_IM_boolean(num_intermediate,Post-Pre,Post,At,bt,V,m0,flag_ILP);

            if MILP.exitflag~=1
                if MILP.exitflag == 10
                    simulation.MILP.success = 0;
                    simulation.MILP.exitflag = MILP.exitflag;
                    simulation.MILP.num_intermediate = num_intermediate;
                    flag = 0;
                    return;
                elseif MILP.exitflag == 3600
                    simulation.MILP.success = 0;
                    simulation.MILP.exitflag = MILP.exitflag;
                    simulation.MILP.num_intermediate = num_intermediate;
                    simulation.MILP.runtime = MILP.runtime;
                    flag = 0; %time limit has been exceeded
                    return;
                else
                    if num_intermediate < nR
                        num_intermediate = num_intermediate + 1;
                    else
                        fprintf('Number of intermediary markings is greater than the number of robots.\n');
                        fprintf('Solution not found!\n');
                        simulation.MILP.success = 0;
                        simulation.MILP.exitflag = MILP.exitflag;
                        simulation.MILP.num_intermediate = num_intermediate;
                        flag = 0;
                        return;
                    end
                    simulation.MILP.runtime_all =  simulation.MILP.runtime_all + MILP.runtime;
                end
            else
                % save the data in simulation structure
                simulation.MILP.success = 1;
                simulation.MILP.runtime = MILP.runtime;
                simulation.MILP.cost = 0;
                simulation.MILP.sol = MILP.sol;
                simulation.MILP.num_intermediate = num_intermediate;
                simulation.MILP.exitflag = MILP.exitflag;
                flag_num = 0;
            end
        end
    end

    for i = 1 : num_intermediate
        simulation.MILP.cost = simulation.MILP.cost + sum(MILP.sol((i-1)*(nplaces+ntrans) + nplaces + 1 : i*(nplaces+ntrans)));
        if flag_ILP && simulation.ILP2.success == 1
            simulation.ILP2.cost = simulation.ILP2.cost + sum(ILP2.sol((i-1)*(nplaces+ntrans) + nplaces + 1 : i*(nplaces+ntrans)));
        end
    end

    if flag_ILP
        if simulation.MILP.cost ~= simulation.ILP2.cost
            fprintf(1,'===============\n');
            fprintf(1,'(8) - MILP (m_i, sigma_i, x are unknowns) cost: %i\n',simulation.MILP.cost);
            fprintf(1,'(8) - ILP (m_i, sigma_i and x are unknowns) cost: %i\n',simulation.ILP2.cost);
            fprintf(1,'===============\n');
        else
            fprintf(1,'===============\n');
            fprintf(1,'Optimal cost is obtained! (8) - ILP cost: %i, (8) - MILP cost %i\n',simulation.ILP2.cost,simulation.MILP.cost);
            fprintf(1,'===============\n');
        end
    end
end

fprintf(1,'\nSolution MILP found in %f [secs].',simulation.MILP.runtime);
fprintf(1,'\nOptimal value MILP = %s\n',num2str(simulation.MILP.cost));
fprintf(1,'\nTotal time for solution LP1+MILP: %f [sec].',simulation.LP1.runtime+simulation.MILP.runtime);

if flag_ILP
    fprintf(1,'\nTotal time for solution ILP1+ILP2 with intermediary markings: %f [sec].',simulation.ILP1.runtime+simulation.ILP2.runtime);
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
            current_marking = round(MILP.sol((nplaces+ntrans)*(i-2)+1:(nplaces+ntrans)*(i-2)+nplaces));
        end
        % figure;
        plot_environment(T.rem_cell(find(current_marking)),T.rem_cell(T.props),T.map2D,env_limit,T.Vert); 
        title(sprintf('Efficient path planning (LP1 + MILP): iteration %d',i));
        hold on;

        if num_intermediate == 1
            current_sigma = round(LP1.sol(nplaces+1:nplaces+ntrans));
        else
            current_sigma = round(MILP.sol((nplaces+ntrans)*(i-1)+nplaces+1:(nplaces+ntrans)*(i-1)+nplaces+ntrans));
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

    if flag_ILP && simulation.ILP1.timelimit
        rob_color = hsv(sum(m0));

        for i = 1 : num_intermediate
            if (i == 1)
                current_marking = m0;
            else
                current_marking = round(ILP2.sol((nplaces+ntrans)*(i-2)+1:(nplaces+ntrans)*(i-2)+nplaces));
            end
            %figure;
            plot_environment(T.rem_cell(find(current_marking)),T.rem_cell(T.props),T.map2D,env_limit,T.Vert);
            title(sprintf('ILP formulation: iteration %d',i));
            hold on;

            if num_intermediate == 1
                current_sigma = round(ILP1.sol(nplaces+1:nplaces+ntrans));
            else
                current_sigma = round(ILP2.sol((nplaces+ntrans)*(i-1)+nplaces+1:(nplaces+ntrans)*(i-1)+nplaces+ntrans));
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