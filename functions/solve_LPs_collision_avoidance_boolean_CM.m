function [simulation,flag] = solve_LPs_collision_avoidance_boolean_CM(Post,Pre,At,bt,m0,T,flag_ILP,plot_animation)

flag = 1;
nplaces = size(Post,1);
ntrans  = size(Post,2);
nR      = sum(m0);

% Construcción de V 

np = size(At, 2);
V  = zeros(np, nplaces);

% goals: vector columna con los índices de los lugares finales (1..nplaces)
goals = double(T.props(:));                 % o el vector que uses
idx = sub2ind([np, nplaces], (1:np)', goals);
V(idx) = 1;

% Constantes (evitar "números mágicos")
EXIT_INFEASIBLE = 10;
EXIT_TIMELIMIT  = 3600;

% ====================== PRIMER PROBLEMA (LP1 + ILP1 opcional) ======================
if flag_ILP
    [LP1,ILP1] = solve_mILP_boolean(Post,Pre,V,At,bt,m0,flag_ILP);

    if ILP1.exitflag ~= 1
        if LP1.exitflag == EXIT_INFEASIBLE && ILP1.exitflag == EXIT_INFEASIBLE
            % Ambos infactibles
            simulation.ILP1.timelimit = 0; simulation.ILP1.runtime = 0;
            simulation.ILP1.cost = 0;      simulation.ILP1.sol = 0;
            simulation.ILP1.cellCapacity = 0; simulation.ILP1.exitflag = ILP1.exitflag;

            simulation.LP1.timelimit = 0;  simulation.LP1.runtime = 0;
            simulation.LP1.cost = 0;       simulation.LP1.sol = 0;
            simulation.LP1.cellCapacity = 0; simulation.LP1.exitflag = LP1.exitflag;

            flag = 0; return;
        elseif ILP1.exitflag == EXIT_TIMELIMIT
            simulation.ILP1.timelimit = 0;
            simulation.ILP1.runtime   = ILP1.runtime;
            simulation.ILP1.exitflag  = ILP1.exitflag;
        else
            fprintf('ILP >>> exitflag ILP = %i.\n', ILP1.exitflag);
            simulation.ILP1.timelimit = 0;
            simulation.ILP1.runtime   = ILP1.runtime;
            simulation.ILP1.cost = 0; simulation.ILP1.sol = 0;
            simulation.ILP1.cellCapacity = 0;
            simulation.ILP1.exitflag     = ILP1.exitflag;
        end
    else
        % ILP1 OK
        simulation.ILP1.timelimit    = 1;
        simulation.ILP1.runtime      = ILP1.runtime;
        simulation.ILP1.cost         = sum(ILP1.sol(nplaces+1:nplaces+ntrans));
        simulation.ILP1.sol          = ILP1.sol;
        simulation.ILP1.cellCapacity = ILP1.sol(end);
        simulation.ILP1.exitflag     = ILP1.exitflag;
    end
else
    LP1 = solve_mILP_boolean(Post,Pre,V,At,bt,m0,flag_ILP);
end

% LP1: gestionar salida
if LP1.exitflag ~= 1
    if LP1.exitflag == EXIT_INFEASIBLE
        simulation.LP1.timelimit = 0; simulation.LP1.runtime = 0;
        simulation.LP1.exitflag  = LP1.exitflag;
        simulation.LP1.cost = 0; simulation.LP1.sol = 0; simulation.LP1.cellCapacity = 0;
    elseif LP1.exitflag == EXIT_TIMELIMIT
        simulation.LP1.timelimit = 0; simulation.LP1.runtime = LP1.runtime;
        simulation.LP1.exitflag  = LP1.exitflag;
    else
        fprintf('LP1 >>> exitflag LP1 = %i.\n', LP1.exitflag);
        simulation.LP1.timelimit = 0; simulation.LP1.runtime = LP1.runtime;
        simulation.LP1.exitflag  = LP1.exitflag;
        simulation.LP1.cost = 0; simulation.LP1.sol = 0; simulation.LP1.cellCapacity = 0;
    end
    flag = 0; return;
else
    % LP1 OK
    x = LP1.sol(nplaces+ntrans+1 : nplaces+ntrans+np);
    s = LP1.sol(end);  % cell capacity

    simulation.LP1.timelimit    = 1;
    simulation.LP1.runtime      = LP1.runtime;
    simulation.LP1.cost         = sum(LP1.sol(nplaces+1:nplaces+ntrans));
    simulation.LP1.sol          = LP1.sol;
    simulation.LP1.cellCapacity = LP1.sol(end);
    simulation.LP1.exitflag     = LP1.exitflag;
end

% Prints MILP1
fprintf(1,'\n============================================================');
fprintf(1,'\n    Solution for MILP1 and ILP1 obtained using intlinprog solver ');
fprintf(1,'\n============================================================');
fprintf(1,'\nSolution for MILP1 found in %f [secs].', simulation.LP1.runtime);
fprintf(1,'\nOptimal solution MILP1 = %s', num2str(simulation.LP1.cost));
fprintf(1,'\nInfinite norm for MILP1 of Post * sigma = %s', num2str(simulation.LP1.cellCapacity));
if flag_ILP
    fprintf(1,'\nSolution for ILP1 found in %f [secs].', simulation.ILP1.runtime);
    fprintf(1,'\nOptimal solution ILP1 = %s', num2str(simulation.ILP1.cost));
    fprintf(1,'\nInfinite norm for ILP1 of Post * sigma = %s', num2str(simulation.ILP1.cellCapacity));
end
% ====================== DECISIÓN DE SEGUNDO PROBLEMA ======================
if norm(s - 1) < 1000*eps && norm(x - round(x)) < 1000*eps
    % Capacidad 1 y x entero -> no resolver segundo problema
    num_intermediate = 1;
    fprintf("\nCollision avoidance is already imposed. The second problem will not be solved.\n");

    % Estructuras triviales
    simulation.MILP.runtime = 0; simulation.MILP.cost = 0;
    simulation.MILP.sol = 0;     simulation.MILP.num_intermediate = 0;
    if flag_ILP
        simulation.ILP2.runtime = 0; simulation.ILP2.cost = 0;
    end

else
    fprintf(1,'\n\n=============================================');
    fprintf(1,'\n           Start solving second problem...');
    fprintf(1,'\n=============================================');

    clear Aeq1; clear beq1; clear Aineq1; clear bineq1;

    num_intermediate = ceil(s);
    if num_intermediate < 1, num_intermediate = 1; end

    flag_num = 1;
    simulation.ILP2.runtime_all = 0;
    simulation.MILP.runtime_all = 0;
    simulation.ILP2.cost = 0;
    simulation.MILP.cost  = 0;

    while flag_num
        fprintf(1,'\nSolve second problem with %d intermediate markings\n', num_intermediate);

        if flag_ILP
            [MILP,ILP2] = solve_LP_IM_boolean(num_intermediate, Post-Pre, Post, At, bt, V, m0, flag_ILP);

            % ---- ILP2 ----
            if ILP2.exitflag ~= 1
                if ILP2.exitflag == EXIT_INFEASIBLE && MILP.exitflag == EXIT_INFEASIBLE
                    simulation.ILP2.success = 0; simulation.MILP.success = 0;
                    simulation.ILP2.exitflag = ILP2.exitflag; simulation.MILP.exitflag = MILP.exitflag;
                    simulation.ILP2.num_intermediate = num_intermediate; simulation.MILP.num_intermediate = num_intermediate;
                    flag = 0; return;

                elseif ILP2.exitflag == EXIT_TIMELIMIT
                    simulation.ILP2.timelimit = EXIT_TIMELIMIT;
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
                        simulation.ILP2.success = 0; simulation.MILP.success = 0;
                        simulation.ILP2.exitflag = ILP2.exitflag; simulation.MILP.exitflag = MILP.exitflag;
                        flag = 0; return;
                    end
                    simulation.ILP2.runtime_all = simulation.ILP2.runtime_all + ILP2.runtime;
                    simulation.MILP.runtime_all = simulation.MILP.runtime_all + MILP.runtime;
                end
            else
                if flag_ILP && exist('ILP1','var') && (ILP1.runtime + ILP2.runtime < EXIT_TIMELIMIT)
                    flag_num = 0;
                    simulation.ILP2.success = 1; simulation.ILP2.runtime = ILP2.runtime;
                    simulation.ILP2.cost = 0;   simulation.ILP2.sol = ILP2.sol;
                    simulation.ILP2.exitflag = ILP2.exitflag; simulation.ILP2.num_intermediate = num_intermediate;
                else
                    flag_num = 0;
                    simulation.ILP2.success = EXIT_TIMELIMIT;
                    simulation.ILP2.runtime = ILP2.runtime;
                    simulation.ILP2.exitflag = ILP2.exitflag;
                    simulation.ILP2.cost = 0; simulation.ILP2.sol = 0;
                    simulation.ILP2.num_intermediate = num_intermediate;
                end
            end

            % ---- MILP ----
            if MILP.exitflag ~= 1
                if MILP.exitflag == EXIT_INFEASIBLE
                    simulation.MILP.success = 0; simulation.MILP.exitflag = MILP.exitflag;
                    simulation.MILP.num_intermediate = num_intermediate;
                    flag = 0; return;
                elseif MILP.exitflag == EXIT_TIMELIMIT
                    simulation.MILP.success = 0; simulation.MILP.exitflag = MILP.exitflag;
                    simulation.MILP.runtime = MILP.runtime; simulation.MILP.num_intermediate = num_intermediate;
                    flag = 0; return;
                end
            else
                simulation.MILP.success = 1; simulation.MILP.runtime = MILP.runtime;
                simulation.MILP.cost = 0;    simulation.MILP.sol = MILP.sol;
                simulation.MILP.num_intermediate = num_intermediate; simulation.MILP.exitflag = MILP.exitflag;
                flag_num = 0;
            end

        else
            % ---- Sólo MILP ----
            MILP = solve_LP_IM_boolean(num_intermediate, Post-Pre, Post, At, bt, V, m0, flag_ILP);

            if MILP.exitflag ~= 1
                if MILP.exitflag == EXIT_INFEASIBLE
                    simulation.MILP.success = 0; simulation.MILP.exitflag = MILP.exitflag;
                    simulation.MILP.num_intermediate = num_intermediate;
                    flag = 0; return;
                elseif MILP.exitflag == EXIT_TIMELIMIT
                    simulation.MILP.success = 0; simulation.MILP.exitflag = MILP.exitflag;
                    simulation.MILP.num_intermediate = num_intermediate; simulation.MILP.runtime = MILP.runtime;
                    flag = 0; return;
                else
                    if num_intermediate < nR
                        num_intermediate = num_intermediate + 1;
                    else
                        fprintf('Number of intermediary markings is greater than the number of robots.\n');
                        fprintf('Solution not found!\n');
                        simulation.MILP.success = 0; simulation.MILP.exitflag = MILP.exitflag;
                        simulation.MILP.num_intermediate = num_intermediate;
                        flag = 0; return;
                    end
                    simulation.MILP.runtime_all = simulation.MILP.runtime_all + MILP.runtime;
                end
            else
                simulation.MILP.success = 1; simulation.MILP.runtime = MILP.runtime;
                simulation.MILP.cost = 0;    simulation.MILP.sol = MILP.sol;
                simulation.MILP.num_intermediate = num_intermediate; simulation.MILP.exitflag = MILP.exitflag;
                flag_num = 0;
            end
        end
    end

    % Acumular costes
    for i = 1:num_intermediate
        base = (i-1)*(nplaces+ntrans);
        simulation.MILP.cost = simulation.MILP.cost + sum(MILP.sol(base + nplaces + 1 : base + nplaces + ntrans));
        if flag_ILP && isfield(simulation,'ILP2') && isfield(simulation.ILP2,'success') && simulation.ILP2.success == 1
            simulation.ILP2.cost = simulation.ILP2.cost + sum(ILP2.sol(base + nplaces + 1 : base + nplaces + ntrans));
        end
    end

    % Comparativa de costes (si aplica)
    if flag_ILP
        if simulation.MILP.cost ~= simulation.ILP2.cost
            fprintf(1,'===============\n');
            fprintf(1,'(8) - MILP (m_i, sigma_i, x are unknowns) cost: %i\n', simulation.MILP.cost);
            fprintf(1,'(8) - ILP (m_i, sigma_i and x are unknowns) cost: %i\n', simulation.ILP2.cost);
            fprintf(1,'===============\n');
        else
            fprintf(1,'===============\n');
            fprintf(1,'Optimal cost is obtained! (8) - ILP cost: %i, (8) - MILP cost %i\n', simulation.ILP2.cost, simulation.MILP.cost);
            fprintf(1,'===============\n');
        end
    end
end

% ====================== RESÚMENES FINALES ======================
%fprintf(1,'\nSolution MILP found in %f [secs].', simulation.MILP.runtime);
%fprintf(1,'\nOptimal value MILP = %s\n', num2str(simulation.MILP.cost));
fprintf(1,'\nTotal time for solution LP1+MILP: %f [sec].', simulation.LP1.runtime + simulation.MILP.runtime);
if flag_ILP
    fprintf(1,'\nTotal time for solution ILP1+ILP2 with intermediary markings: %f [sec].', simulation.ILP1.runtime + simulation.ILP2.runtime);
end

% ====================== PLOT ======================
if plot_animation
    rob_color = hsv(sum(m0));

    for i = 1:num_intermediate
        if i == 1
            current_marking = m0;
        else
            current_marking = round(simulation.MILP.sol((nplaces+ntrans)*(i-2)+1 : (nplaces+ntrans)*(i-2)+nplaces));
        end

        start_idx = T.rem_cells(find(current_marking));
        goal_idx  = T.rem_cells(T.props);

        [H, W] = size(T.map2D);
        [start_r, start_c] = ind2sub([H, W], start_idx);
        [goal_r,  goal_c ] = ind2sub([H, W], goal_idx);

        selectedStart = [start_c - 0.5 , start_r - 0.5 ];
        selectedFin   = [goal_c  - 0.5 , goal_r  - 0.5 ];

        plot_environment_new_SG_boolean(selectedStart, selectedFin, T.map2D);
        title(sprintf('Efficient path planning (LP1 + LP2): iteration %d', i));
        hold on;

        if num_intermediate == 1
            current_sigma = round(simulation.LP1.sol(nplaces+1 : nplaces+ntrans));
        else
            current_sigma = round(simulation.MILP.sol((nplaces+ntrans)*(i-1)+nplaces+1 : (nplaces+ntrans)*(i-1)+nplaces+ntrans));
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
        %saveas(gcf, sprintf('smart_plant_eff_it%d.fig', i));
    end
end
