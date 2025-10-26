function [simulation,flag] = solve_LPs_collision_avoidance_CM(Post,Pre,mf,m0,T,flag_ILP,plot_animation)
% Near-original implementation, simpler and with fewer repetitions.

flag = 1;
nplaces = size(Post,1);
ntrans  = size(Post,2);
nR      = sum(m0);

EXIT_INFEASIBLE = 10;
EXIT_TIMELIMIT  = 3600;

fprintf(1,'\nSolve first problem\n');

% -------------------- LP1 (+ ILP1 opcional) -------------------- %
if flag_ILP
    [LP1,ILP1] = solve_mILPr(Post,Pre,mf,m0,flag_ILP);

    if ILP1.exitflag ~= 1
        if LP1.exitflag == EXIT_INFEASIBLE && ILP1.exitflag == EXIT_INFEASIBLE
            % Ambos infeasible
            simulation.ILP1.timelimit     = 0;  simulation.ILP1.runtime = 0;
            simulation.ILP1.cost          = 0;  simulation.ILP1.sol     = 0;
            simulation.ILP1.cellCapacity  = 0;  simulation.ILP1.exitflag = ILP1.exitflag;

            simulation.LP1.timelimit      = 0;  simulation.LP1.runtime = 0;
            simulation.LP1.cost           = 0;  simulation.LP1.sol     = 0;
            simulation.LP1.cellCapacity   = 0;  simulation.LP1.exitflag = LP1.exitflag;

            flag = 0; return;
        elseif ILP1.exitflag == EXIT_TIMELIMIT
            % Timeout en ILP1, seguimos con LP1
            simulation.ILP1.timelimit = 0;
            simulation.ILP1.runtime   = ILP1.runtime;
            simulation.ILP1.exitflag  = ILP1.exitflag;
        else
            % Otro fallo en ILP1
            fprintf('ILP >>> exitflag ILP = %i.\n', ILP1.exitflag);
            simulation.ILP1.timelimit     = 0;
            simulation.ILP1.runtime       = ILP1.runtime;
            simulation.ILP1.cost          = 0;
            simulation.ILP1.sol           = 0;
            simulation.ILP1.cellCapacity  = 0;
            simulation.ILP1.exitflag      = ILP1.exitflag;
        end
    else
        % ILP1 OK
        simulation.ILP1.timelimit     = 1;
        simulation.ILP1.runtime       = ILP1.runtime;
        simulation.ILP1.cost          = sum(ILP1.sol(1:ntrans));
        simulation.ILP1.sol           = ILP1.sol(1:ntrans);
        simulation.ILP1.cellCapacity  = ILP1.sol(end);
        simulation.ILP1.exitflag      = ILP1.exitflag;
    end
else
    LP1 = solve_mILPr(Post,Pre,mf,m0,flag_ILP);
end

% LP1: gestionar salida
if LP1.exitflag ~= 1
    if LP1.exitflag == EXIT_INFEASIBLE
        simulation.LP1.timelimit     = 0;  simulation.LP1.runtime = 0;
        simulation.LP1.exitflag      = LP1.exitflag;
        simulation.LP1.cost          = 0;  simulation.LP1.sol     = 0;
        simulation.LP1.cellCapacity  = 0;
    elseif LP1.exitflag == EXIT_TIMELIMIT
        simulation.LP1.timelimit     = 0;
        simulation.LP1.runtime       = LP1.runtime;
        simulation.LP1.exitflag      = LP1.exitflag;
    else
        fprintf('LP1 >>> exitflag LP1 = %i.\n', LP1.exitflag);
        simulation.LP1.timelimit     = 0;
        simulation.LP1.runtime       = LP1.runtime;
        simulation.LP1.exitflag      = LP1.exitflag;
        simulation.LP1.cost          = 0;  simulation.LP1.sol     = 0;
        simulation.LP1.cellCapacity  = 0;
    end
    flag = 0; return;
else
    % LP1 OK
    simulation.LP1.timelimit     = 1;
    simulation.LP1.runtime       = LP1.runtime;
    simulation.LP1.cost          = sum(LP1.sol(1:ntrans));
    simulation.LP1.sol           = LP1.sol(1:ntrans);
    simulation.LP1.cellCapacity  = LP1.sol(end);
    simulation.LP1.exitflag      = LP1.exitflag;

    s = ceil(LP1.sol(end)); % capacidad de celda
end

% -------------------- Prints LP1 / ILP1 -------------------- %
fprintf(1,'\n------------------------------------------------------------');
fprintf(1,'\n    Solution for LP1 obtained using intlinprog solver ');
fprintf(1,'\n------------------------------------------------------------');
fprintf(1,'\nSolution for LP1 found in %f [secs].', simulation.LP1.runtime);
fprintf(1,'\nOptimal solution for LP1 = %s', num2str(simulation.LP1.cost));
fprintf(1,'\nInfinite norm for LP1 of Post * sigma = %s', num2str(simulation.LP1.cellCapacity));
if flag_ILP
    fprintf(1,'\n----------------------------------------------------------');
    fprintf(1,'\nSolution for IP1 found in %f [secs].', simulation.ILP1.runtime);
    fprintf(1,'\nOptimal solution for IP1 = %s', num2str(sum(simulation.ILP1.sol(1:ntrans))));
    fprintf(1,'\nInfinite norm for IP1 of Post * sigma = %s', num2str(simulation.ILP1.cellCapacity));
end

% -------------------- Decisión sobre segundo problema -------------------- %
if norm(s - 1) < 1000*eps
    fprintf('\nCell capacity is equal with 1. Collision-free solution obtained!\n');
else
    fprintf(1,'\n\n=============================================');
    fprintf(1,'\n           Start solving second problem...');
    fprintf(1,'\n=============================================');
end

clear Aeq1; clear beq1; clear Aineq1; clear bineq1;
num_intermediate = s;

% -------------------- LP2 (+ ILP2 opcional) -------------------- %
if num_intermediate > 1
    flag_num = 1;
    simulation.ILP2.runtime_all = 0;
    simulation.LP2.runtime_all  = 0;
    simulation.ILP2.cost = 0;
    simulation.LP2.cost  = 0;

    while flag_num
        fprintf(1,'\nSolve second problem with %d intermediate markings\n', num_intermediate);

        if flag_ILP
            [LP2,ILP2] = solve_LP_IMr(num_intermediate, Post-Pre, Post, m0, mf, flag_ILP);

            % ---- ILP2 ----
            if ILP2.exitflag ~= 1
                if ILP2.exitflag == EXIT_INFEASIBLE && LP2.exitflag == EXIT_INFEASIBLE
                    simulation.ILP2.success = 0;   simulation.LP2.success = 0;
                    simulation.ILP2.exitflag = ILP2.exitflag;
                    simulation.LP2.exitflag  = LP2.exitflag;
                    simulation.ILP2.num_intermediate = num_intermediate;
                    simulation.LP2.num_intermediate  = num_intermediate;
                    flag = 0; return;

                elseif ILP2.exitflag == EXIT_TIMELIMIT
                    simulation.ILP2.success  = 0;
                    simulation.ILP2.runtime  = ILP2.runtime;
                    simulation.ILP2.exitflag = ILP2.exitflag;
                    simulation.ILP2.num_intermediate = num_intermediate;

                else
                    if num_intermediate < nR
                        num_intermediate = num_intermediate + 1;
                    else
                        simulation.ILP2.num_intermediate = num_intermediate;
                        simulation.LP2.num_intermediate  = num_intermediate;
                        simulation.ILP2.success = 0;   simulation.LP2.success = 0;
                        simulation.ILP2.exitflag = ILP2.exitflag;
                        simulation.LP2.exitflag  = LP2.exitflag;
                        flag = 0; return;
                    end
                    simulation.ILP2.runtime_all = simulation.ILP2.runtime_all + ILP2.runtime;
                    simulation.LP2.runtime_all  = simulation.LP2.runtime_all  + LP2.runtime;
                end
            else
                if ILP1.runtime + ILP2.runtime < EXIT_TIMELIMIT
                    flag_num = 0;
                    simulation.ILP2.success  = 1;
                    simulation.ILP2.runtime  = ILP2.runtime;
                    simulation.ILP2.cost     = 0;
                    simulation.ILP2.sol      = ILP2.sol;
                    simulation.ILP2.exitflag = ILP2.exitflag;
                    simulation.ILP2.num_intermediate = num_intermediate;
                else
                    flag_num = 0;
                    simulation.ILP2.success  = EXIT_TIMELIMIT;
                    simulation.ILP2.runtime  = ILP2.runtime;
                    simulation.ILP2.exitflag = ILP2.exitflag;
                    simulation.ILP2.cost     = 0;
                    simulation.ILP2.sol      = 0;
                    simulation.ILP2.num_intermediate = num_intermediate;
                end
            end

            % ---- LP2 ----
            if LP2.exitflag ~= 1
                if LP2.exitflag == EXIT_INFEASIBLE
                    simulation.LP2.success = 0;
                    simulation.LP2.exitflag = LP2.exitflag;
                    simulation.LP2.num_intermediate = num_intermediate;
                    flag = 0; return;

                elseif LP2.exitflag == EXIT_TIMELIMIT
                    simulation.LP2.exitflag = LP2.exitflag;
                    simulation.LP2.runtime  = LP2.runtime;
                    simulation.LP2.num_intermediate = num_intermediate;
                    simulation.LP2.success  = EXIT_TIMELIMIT;
                    flag = 0; return;
                end
            else
                simulation.LP2.success = 1;
                simulation.LP2.runtime = LP2.runtime;
                simulation.LP2.cost    = 0;
                simulation.LP2.sol     = LP2.sol;
                simulation.LP2.no_inter_markings = num_intermediate;
                simulation.LP2.exitflag = LP2.exitflag;
                flag_num = 0;
            end

        else
            % ---- Solo LP2 ----
            LP2 = solve_LP_IMr(num_intermediate, Post-Pre, Post, m0, mf, flag_ILP);

            if LP2.exitflag ~= 1
                if LP2.exitflag == EXIT_INFEASIBLE
                    simulation.LP2.success = 0;
                    simulation.LP2.exitflag = LP2.exitflag;
                    simulation.LP2.num_intermediate = num_intermediate;
                    flag = 0; return;

                elseif LP2.exitflag == EXIT_TIMELIMIT
                    simulation.LP2.exitflag = LP2.exitflag;
                    simulation.LP2.runtime  = LP2.runtime;
                    simulation.LP2.num_intermediate = num_intermediate;
                    simulation.LP2.success  = 0;
                    flag = 0; return;

                else
                    if num_intermediate < nR
                        num_intermediate = num_intermediate + 1;
                    else
                        fprintf('Number of intermediary markings is greater than the number of robots.\n');
                        fprintf('Solution not found!\n');
                        simulation.LP2.num_intermediate = num_intermediate;
                        simulation.LP2.success = 0;
                        simulation.LP2.exitflag = LP2.exitflag;
                        flag = 0; return;
                    end
                    simulation.LP2.runtime_all = simulation.LP2.runtime_all + LP2.runtime;
                end
            else
                simulation.LP2.success = 1;
                simulation.LP2.runtime = LP2.runtime;
                simulation.LP2.cost    = 0;
                simulation.LP2.sol     = LP2.sol;
                simulation.LP2.no_inter_markings = num_intermediate;
                simulation.LP2.exitflag = LP2.exitflag;
                flag_num = 0;
            end
        end
    end

    % Acumular costes LP2/ILP2
    for i = 1:num_intermediate
        baseIdx = (i-1)*(nplaces+ntrans);
        simulation.LP2.cost = simulation.LP2.cost + sum(LP2.sol(baseIdx + nplaces + 1 : baseIdx + nplaces + ntrans));
        if flag_ILP && isfield(simulation,'ILP2') && isfield(simulation.ILP2,'success') && simulation.ILP2.success == 1
            simulation.ILP2.cost = simulation.ILP2.cost + sum(ILP2.sol(baseIdx + nplaces + 1 : baseIdx + nplaces + ntrans));
        end
    end
else
    % No se resuelve segundo problema (trivial)
    simulation.LP2.runtime = 0;
    simulation.LP2.cost    = 0;
    simulation.LP2.sol     = 0;
    simulation.LP2.no_inter_markings = 0;

    if flag_ILP
        simulation.ILP2.runtime = 0;
        simulation.ILP2.cost    = 0;
    end
end

% -------------------- Prints finales -------------------- %
if num_intermediate > 1
    fprintf(1,'\nSolution LP2 found in %f [secs].', simulation.LP2.runtime);
    fprintf(1,'\nOptimal value LP2 = %s', num2str(simulation.LP2.cost));
    fprintf(1,'\nTotal time for solution LP1+LP2: %f [sec].', simulation.LP1.runtime + simulation.LP2.runtime);
    if flag_ILP
        fprintf(1,'\n\nSolution IP2 found in %f [secs].', simulation.ILP2.runtime);
        fprintf(1,'\nOptimal value IP2 = %s', num2str(simulation.ILP2.cost));
        fprintf(1,'\nTotal time for solution ILP1+ILP2 with intermediary markings: %f [sec].', simulation.ILP1.runtime + simulation.ILP2.runtime);
    end
else
    fprintf(1,'\nTotal time for solution LP1+LP2: %f [sec].', simulation.LP1.runtime + simulation.LP2.runtime);
    fprintf(1,'\nTotal time for solution ILP1+ILP2 with intermediary markings: %f [sec].', simulation.ILP1.runtime + simulation.ILP2.runtime);
end

% -------------------- Plot -------------------- %
if plot_animation
    rob_color = hsv(sum(m0));

    for i = 1:num_intermediate
        if i == 1
            current_marking = m0;
        else
            current_marking = round(LP2.sol((nplaces+ntrans)*(i-2)+1:(nplaces+ntrans)*(i-2)+nplaces));
        end

        start_idx = T.rem_cells(find(current_marking));
        goal_idx  = T.rem_cells(T.props);

        n_r = min(numel(start_idx), numel(goal_idx));
        if n_r == 0
            warning('No start/goal pairs to plot.');
            return;
        end
        start_idx = start_idx(1:n_r);
        goal_idx  = goal_idx(1:n_r);

        [H, W] = size(T.map2D);
        [start_r, start_c] = ind2sub([H, W], start_idx(1:n_r));
        [goal_r,  goal_c ] = ind2sub([H, W], goal_idx (1:n_r));

        selectedStart = [start_c - 0.5 , start_r - 0.5 ];
        selectedFin   = [goal_c  - 0.5 , goal_r  - 0.5 ];

        plot_environment_new_SG(selectedStart, selectedFin, T.map2D, T);
        title(sprintf('Efficient path planning (LP1 + LP2): iteration %d', i));
        hold on;

        if num_intermediate == 1
            current_sigma = round(LP1.sol(1:ntrans));
        else
            base = (nplaces+ntrans)*(i-1);
            current_sigma = round(LP2.sol(base + nplaces + 1 : base + nplaces + ntrans));
        end

        try
            [feasible_sigma, Rob_places, ~, ~] = sigma2trajectories(Pre,Post,current_marking,current_sigma,find(current_marking));
            if ~feasible_sigma
                error('something wrong happened!');
            end
        catch
            warning('something wrong happened!');
            flag = 0; return;
        end

        nRloc = numel(Rob_places);
        if nRloc == 0, continue; end
        if size(rob_color,1) < nRloc, rob_color = hsv(nRloc); end

        for j = 1:nRloc
            traj = Rob_places{j};
            XY = places2xy(traj, T);
            plot(XY(:,1), XY(:,2), '-', 'LineWidth', 1.5, 'Color', rob_color(j,:));
            c_start = XY(1,:); c_goal = XY(end,:);
            scatter(c_start(1), c_start(2), 36, rob_color(j,:), '^', 'filled', 'MarkerEdgeColor', rob_color(j,:));
            scatter(c_goal(1),  c_goal(2),  36, rob_color(j,:), 'd', 'filled', 'MarkerEdgeColor',  rob_color(j,:));
        end
    end
end
