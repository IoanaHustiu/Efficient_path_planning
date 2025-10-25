function [simulation,flag] = solve_LPs_collision_avoidance(Post,Pre,mf,m0,T,flag_ILP,plot_animation)
flag = 1;
nplaces = size(Post,1);
ntrans = size(Post,2);
nR = sum(m0);
np = nR;
env_limit = size(T.map2D);

fprintf(1,'\nSolve first problem\n');
if flag_ILP
    [LP1,ILP1] = solve_mILPr(Post,Pre,mf,m0,flag_ILP);

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
        simulation.ILP1.cost = sum(ILP1.sol(1:ntrans));
        simulation.ILP1.sol = ILP1.sol(1:ntrans);
        simulation.ILP1.cellCapacity = ILP1.sol(end);
        simulation.ILP1.exitflag = ILP1.exitflag;
        s_ILP = ceil(ILP1.sol(end)); % cell capacity
    end

else
    LP1 = solve_mILPr(Post,Pre,mf,m0,flag_ILP);
end

if LP1.exitflag~=1
    if LP1.exitflag == 10
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
    % save the data in simulation structure
    simulation.LP1.timelimit = 1;
    simulation.LP1.runtime = LP1.runtime;
    simulation.LP1.cost = sum(LP1.sol(1:ntrans));
    simulation.LP1.sol = LP1.sol(1:ntrans);
    simulation.LP1.cellCapacity = LP1.sol(end);
    simulation.LP1.exitflag = LP1.exitflag;

    s = ceil(LP1.sol(end)); % cell capacity
end

%check if the cell capacity is equal with 1
if norm(s - 1) < 1000*eps %cell capacity is 1 and collision avoidance is not required
    fprintf('Cell capacity is equal with 1. Collision-free solution obtained!\n');

    fprintf(1,'\n============================================================');
    fprintf(1,'\n    Solution for LP1 obtained using intlinprog solver ');
    fprintf(1,'\n============================================================');
    fprintf(1,'\nSolution for LP1 found in %f [secs].',LP1.runtime);
    fprintf(1,'\nOptimal solution for LP1 = %s',num2str(LP1.cost));
    fprintf(1,'\nInfinite norm for LP1 of Post * sigma = %s',num2str(LP1.sol(end)));
end

if plot_animation
    fprintf(1,'\n============================================================');
    fprintf(1,'\n    Solution for LP1 obtained using intlinprog solver ');
    fprintf(1,'\n============================================================');
    fprintf(1,'\nSolution for LP1 found in %f [secs].',LP1.runtime);
    fprintf(1,'\nOptimal solution for LP1 = %s',num2str(LP1.cost));
    fprintf(1,'\nInfinite norm for LP1 of Post * sigma = %s',num2str(LP1.sol(end)));
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
    simulation.ILP2.runtime_all = 0;
    simulation.LP2.runtime_all = 0;
    simulation.ILP2.cost = 0;
    simulation.LP2.cost = 0;

    while flag_num
        fprintf(1,'\nSolve second problem with %d intermediate markings\n',num_intermediate);

        if flag_ILP
            [LP2,ILP2] = solve_LP_IMr(num_intermediate,Post-Pre,Post,m0,mf,flag_ILP);

            if ILP2.exitflag~=1
                if ILP2.exitflag == 10 && LP2.exitflag == 10
                    simulation.ILP2.success = 0;
                    simulation.LP2.success = 0;
                    simulation.ILP2.exitflag = ILP2.exitflag;
                    simulation.LP2.exitflag = LP2.exitflag;
                    simulation.ILP2.num_intermediate = num_intermediate;
                    simulation.LP2.num_intermediate = num_intermediate;

                    flag = 0;
                    flag_num = 0;
                    return;
                elseif ILP2.exitflag == 3600
                    simulation.ILP2.success = 0;
                    simulation.ILP2.runtime = ILP2.runtime;
                    simulation.ILP2.exitflag = ILP2.exitflag;
                    simulation.ILP2.num_intermediate = num_intermediate;
                else
                    if num_intermediate < nR
                        num_intermediate = num_intermediate + 1;
                    else
                        simulation.ILP2.num_intermediate = num_intermediate;
                        simulation.LP2.num_intermediate = num_intermediate;
                        simulation.ILP2.success = 0;
                        simulation.LP2.success = 0;
                        simulation.ILP2.exitflag = ILP2.exitflag;
                        simulation.LP2.exitflag = LP2.exitflag;

                        flag = 0;
                        flag_num = 0;
                        return;
                    end
                    simulation.ILP2.runtime_all = simulation.ILP2.runtime_all + ILP2.runtime;
                    simulation.LP2.runtime_all =  simulation.LP2.runtime_all + LP2.runtime;
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
                    % flag = 0;
                    % return;
                end
            end
                if LP2.exitflag~=1
                    if LP2.exitflag == 10
                        simulation.LP2.success = 0;
                        simulation.LP2.exitflag = LP2.exitflag;
                        simulation.LP2.num_intermediate = num_intermediate;
                        flag = 0;
                        return;
                    elseif LP2.exitflag == 3600
                        simulation.LP2.exitflag = LP2.exitflag;
                        simulation.LP2.runtime = LP2.runtime;
                        simulation.LP2.num_intermediate = num_intermediate;
                        simulation.LP2.success = 3600;
                        flag = 0; %time limit has been exceeded
                        return;
                    end
                else
                    % save the data in simulation structure
                    simulation.LP2.success = 1;
                    simulation.LP2.runtime = LP2.runtime;
                    simulation.LP2.cost = 0;
                    simulation.LP2.sol = LP2.sol;
                    simulation.LP2.no_inter_markings = num_intermediate;
                    simulation.LP2.exitflag = LP2.exitflag;
                    flag_num = 0;
                end
            % end
        else
            LP2 = solve_LP_IMr(num_intermediate,Post-Pre,Post,m0,mf,flag_ILP);

            if LP2.exitflag~=1
                if LP2.exitflag == 10
                    simulation.LP2.success = 0;
                    simulation.LP2.exitflag = LP2.exitflag;
                    simulation.LP2.num_intermediate = num_intermediate;
                    flag = 0;
                    return;
                elseif LP2.exitflag == 3600
                    simulation.LP2.exitflag = LP2.exitflag;
                    simulation.LP2.runtime = LP2.runtime;
                    simulation.LP2.num_intermediate = num_intermediate;
                    simulation.LP2.success = 0;
                    flag = 0; %time limit has been exceeded
                    return;
                else
                    if num_intermediate < nR
                        num_intermediate = num_intermediate + 1;
                    else
                        fprintf('Number of intermediary markings is greater than the number of robots.\n');
                        fprintf('Solution not found!\n');
                        simulation.LP2.num_intermediate = num_intermediate;
                        simulation.LP2.success = 0;
                        simulation.LP2.exitflag = LP2.exitflag;
                        flag = 0;
                        return;
                    end
                    simulation.LP2.runtime_all =  simulation.LP2.runtime_all + LP2.runtime;
                end
            else
                % save the data in simulation structure
                simulation.LP2.success = 1;
                simulation.LP2.runtime = LP2.runtime;
                simulation.LP2.cost = 0;
                simulation.LP2.sol = LP2.sol;
                simulation.LP2.no_inter_markings = num_intermediate;
                simulation.LP2.exitflag = LP2.exitflag;
                flag_num = 0;
            end
        end
    end

    for i = 1 : num_intermediate
        simulation.LP2.cost = simulation.LP2.cost + sum(LP2.sol((i-1)*(nplaces+ntrans) + 1 + nplaces:i*(nplaces+ntrans)));
        if flag_ILP && simulation.ILP2.success == 1
            simulation.ILP2.cost = simulation.ILP2.cost + sum(ILP2.sol((i-1)*(nplaces+ntrans) + 1 + nplaces:i*(nplaces+ntrans)));
        end
    end
else
    fprintf("\nCollision avoidance is already imposed. The second problem will not be solved.\n");
    % save the data in simulation structure
    simulation.LP2.runtime = 0;
    simulation.LP2.cost = 0;
    simulation.LP2.sol = 0;
    simulation.LP2.no_inter_markings = 0;

    if flag_ILP
        % save the data in simulation structure
        simulation.ILP2.runtime = 0;
        simulation.ILP2.cost = 0;
    end
end

fprintf(1,'\nSolution LP2 found in %f [secs].',simulation.LP2.runtime);
fprintf(1,'\nOptimal value LP2 = %s\n',num2str(simulation.LP2.cost));
fprintf(1,'\nTotal time for solution LP1+LP2: %f [sec].',simulation.LP1.runtime+simulation.LP2.runtime);

if flag_ILP
    fprintf(1,'\nTotal time for solution ILP1+ILP2 with intermediary markings: %f [sec].',simulation.ILP1.runtime+simulation.ILP2.runtime);
end

if plot_animation
    %rob_color={'b','g','c','m','k','y'};
    % rob_color = hsv(sum(m0));
    rob_color = {'k'};

    for i = 1 : num_intermediate
        if (i == 1)
            current_marking = m0;
        else
            current_marking = round(LP2.sol((nplaces+ntrans)*(i-2)+1:(nplaces+ntrans)*(i-2)+nplaces));
        end
        %figure;
        % plot_environment(T.rem_cell(find(current_marking)),T.rem_cell(T.props),T.map2D,env_limit,T.Vert);
        plot_environment_no_colors(T.rem_cell(find(current_marking)),T.rem_cell(T.props),T.map2D,env_limit,T.Vert); 
        title(sprintf('Efficient path planning (LP1 + LP2): iteration %d',i));
        hold on;

        if num_intermediate == 1
            current_sigma = round(LP1.sol(1:ntrans));
        else
            current_sigma = round(LP2.sol((nplaces+ntrans)*(i-1) + nplaces + 1 : (nplaces+ntrans)*(i-1) + nplaces + ntrans));
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
                    [T.centr{T.rem_cell(traj(k))}(2) T.centr{T.rem_cell(traj(k+1))}(2)],'-','LineWidth',1,'Color','k');
            end
            if (length(traj) > 1)
               init = T.centr{T.rem_cell(traj(1))};
               fill([init(1)-0.15 init(1)+0.15 init(1)+0.15 init(1)-0.15],[init(2)-0.15 init(2)-0.15 init(2)+0.15 init(2)+0.15],'k','EdgeColor','none');
               init = T.centr{T.rem_cell(traj(end))};
               fill([init(1)-0.15 init(1)+0.15 init(1)+0.15 init(1)-0.15],[init(2)-0.15 init(2)-0.15 init(2)+0.15 init(2)+0.15],'k','EdgeColor','none');
            end
            %final destinations already reached
            intersectt = intersect(traj(end),find(mf));
            if ~isempty(intersectt)
                reached = T.centr{T.rem_cell(intersectt)};
                %text(reached(1), reached(2), 'X', 'HorizontalAlignment','center','FontSize',12);
                c = reached;
                theta = linspace(0,2*pi,64);
                r = 0.5;
                x = c(1) + r*cos(theta);
                y = c(2) + r*sin(theta);
                plot(x, y, 'g-', 'LineWidth', 1.5);
            end
        end
    end
    % hold on;
    % plot([-0.55 10.55 10.55 -0.55 -0.55],[-0.55 -0.55 10.55 10.55 -0.55],'k','LineWidth',0.05);
    % axis off;
    % box on;
    % ax1 = gcf;
    % exportgraphics(ax1,'TAPF_it2.png','Resolution',300);
               
    if flag_ILP
        % rob_color = hsv(sum(m0));
        rob_color = {'k'};

        for i = 1 : num_intermediate
            if (i == 1)
                current_marking = m0;
            else
                current_marking = round(ILP2.sol((nplaces+ntrans)*(i-2)+1:(nplaces+ntrans)*(i-2)+nplaces));
            end
            %figure;
            % plot_environment(T.rem_cell(find(current_marking)),T.rem_cell(T.props),T.map2D,env_limit,T.Vert);
            plot_environment_no_colors(T.rem_cell(find(current_marking)),T.rem_cell(T.props),T.map2D,env_limit,T.Vert);
            title(sprintf('ILP formulation: iteration %d',i));
            hold on;

            if num_intermediate == 1
                current_sigma = round(ILP1.sol(1:ntrans));
            else
                current_sigma = round(ILP2.sol((nplaces+ntrans)*(i-1) + nplaces + 1 : (nplaces+ntrans)*(i-1) + nplaces + ntrans));
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
                        [T.centr{T.rem_cell(traj(k))}(2) T.centr{T.rem_cell(traj(k+1))}(2)],'-','LineWidth',1,'Color','k');
                end
                if (length(traj) > 1)
                    init = T.centr{T.rem_cell(traj(1))};
                    fill([init(1)-0.15 init(1)+0.15 init(1)+0.15 init(1)-0.15],[init(2)-0.15 init(2)-0.15 init(2)+0.15 init(2)+0.15],'k','EdgeColor','none');
                    init = T.centr{T.rem_cell(traj(end))};
                    fill([init(1)-0.15 init(1)+0.15 init(1)+0.15 init(1)-0.15],[init(2)-0.15 init(2)-0.15 init(2)+0.15 init(2)+0.15],'k','EdgeColor','none');
                end
                %final destinations already reached
                intersectt = intersect(traj(end),find(mf));
                if ~isempty(intersectt)
                    reached = T.centr{T.rem_cell(intersectt)};
                    %text(reached(1), reached(2), 'X', 'HorizontalAlignment','center','FontSize',12);
                    c = reached;
                    theta = linspace(0,2*pi,64);
                    r = 0.5;
                    x = c(1) + r*cos(theta);
                    y = c(2) + r*sin(theta);
                    plot(x, y, 'g-', 'LineWidth', 1);
                end
            end
        end
    end

end