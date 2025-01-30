function [simulation,flag] = solve_LPs(Post,Pre,m0,mf,T,env_limit,plot_animation,flag_ILP)

% for each simulation the runtime, the cost, the number of iteration/intermediary markings and the cell capacity will be saved
simulation.completeMILP = [];
simulation.completeLP = [];
simulation.interMILP = [];
simulation.interLP = [];
flag = 1;
C = Post - Pre;
nplaces = size(Post,1);
ntrans = size(Post,2);
if (plot_animation)
    plot_environment(T.rem_cell(find(m0)),T.rem_cell(find(mf)),T.obstacles,env_limit,T.Vert); %partition
    title(sprintf('Iteration %d',1));
end

fprintf(1,'Solve first problem\n');

%%%%first problem
f1 = [ones(1,ntrans) sum(m0)*100];
%state equation
Aeq1 = [C zeros(nplaces,1)];
beq1 = mf-m0;
Aineq1 = [Post -ones(nplaces,1)];
bineq1 = -m0;

opt = optimoptions('intlinprog','Display','none','MaxTime',3600);

tic;
[sol1,fval1,exitflag1] = intlinprog(f1,[],Aineq1,bineq1,Aeq1,beq1,zeros(1,size(Aeq1,2)),[],[],opt);
time_LP1 = toc;
if (exitflag1 ~= 1)
    fprintf('Error solving first LP problem!\n');
    fprintf('Generate new initial anf final marking\n');
    flag = 0;
    return;
else
    fprintf('First LP problem solved!\n');    
end

if flag_ILP
    tic;
    [sol1i,fval1i,exitflag1] = intlinprog(f1,1:size(Aeq1,2),Aineq1,bineq1,Aeq1,beq1,zeros(1,size(Aeq1,2)),[],[],opt);
    time_IP1 = toc;

    if (exitflag1 ~= 1)
        fprintf('Error solving first IP problem!\n');
        return;
    else
        fprintf('First IP problem solved!\n');
    end
end


% save the data in simulation structure
if flag_ILP
    simulation.completeMILP.runtime = time_IP1;
    simulation.completeMILP.cost = sum(sol1i(1:ntrans));
    simulation.completeMILP.sol = sol1i(1:end-1);
    simulation.completeMILP.cellCapacity = sol1i(end);
end
simulation.completeLP.runtime = time_LP1;
simulation.completeLP.cost = sum(sol1(1:ntrans));
simulation.completeLP.sol = sol1(1:end-1);
simulation.completeLP.cellCapacity = sol1(end);
simulation.completeLPfixingB.runtime = [];
simulation.completeLPfixingB.cost = [];
simulation.completeLPfixingB.cellCapacity = [];


%check if the solution of the first LP is integer
if norm(ceil(sol1(end)-eps*1000)- sol1(end)) > 1000*eps %solution is not integer
    %fix s
    fprintf(1,'Solve first problem with ceil(s)\n');
    Aeq1 = [Aeq1; [zeros(1,ntrans) 1]];
    beq1 = [beq1; ceil(sol1(end)-eps*1000)];
    tic;
    [sol1e,~,exitflag1] = intlinprog(f1,[],Aineq1,bineq1,Aeq1,beq1,zeros(1,size(Aeq1,2)),[],[],opt);
    time_LP1plus = toc;
    if (exitflag1 ~= 1)
        fprintf('Error solving LP problem with ceil!\n');
        return;
    else
        fprintf('First LP problem with ceil solved!\n');
    end
    simulation.completeLPfixingB.runtime = time_LP1plus;
    simulation.completeLPfixingB.cost = sum(sol1e(1:ntrans));
    simulation.completeLPfixingB.cellCapacity = ceil(sol1(end)-eps*1000);
end


if plot_animation
    fprintf(1,'\n============================================================');
    fprintf(1,'\n    Solution for LP obtained using intlinprog solver ');
    fprintf(1,'\n============================================================');
    fprintf(1,'\nSolution for the first LP found in %f [secs].',time_LP1);
    fprintf(1,'\nOptimal solution for the first LP = %s',num2str(fval1));
    fprintf(1,'\nInfinite norm for the first LP of Post * sigma=%s',num2str(sol1(end)));

    if flag_ILP
        fprintf(1,'\n\nSolution IPL found in %f [secs].',time_IP1);
        fprintf(1,'\nOptimal solution IPL = %s',num2str(fval1i));
        fprintf(1,'\nInfinite norm IPL of Post * sigma=%s',num2str(sol1i(end)));
    end
    %solve the problem with intermediate markings
    fprintf(1,'\n\n=============================================');
    fprintf(1,'\n           Start solving second problem...');
    fprintf(1,'\n=============================================');
end

clear Aeq1;clear beq1;
clear Aineq1;clear bineq1;

exitflag2 = 0;
num_intermediate = ceil(sol1(end)-eps*1000);

if (num_intermediate > 1)
    fprintf(1,'\nSolve second problem with %d intermediate markings \n',num_intermediate);
    while (exitflag2 ~= 1)
        Aeq2 = zeros(num_intermediate*nplaces+1,num_intermediate*(nplaces+ntrans));
        Aineq2 = zeros(num_intermediate*nplaces,num_intermediate*(nplaces+ntrans));
        if (num_intermediate > 1) % b* > 1
            f2 = [zeros(1,nplaces) ones(1,ntrans)];
            Aeq2(1:nplaces,1:nplaces+ntrans) = [eye(nplaces) -C];
            beq2= m0;
            Aineq2(1:nplaces,1:nplaces+ntrans) = [zeros(nplaces,nplaces) Post];
            bineq2 = ones(nplaces,1)-m0;
            bloque=[-eye(nplaces) zeros(nplaces,ntrans) eye(nplaces) -C];
            bloque2=[eye(nplaces) zeros(nplaces,ntrans+nplaces) Post];
            for i = 2 : num_intermediate
                f2 = [f2 zeros(1,nplaces) (num_intermediate - i + 2)*10*sum(m0)*ones(1,ntrans)];
                % add state equation
                Aeq2((i-1)*(nplaces)+1: i*nplaces, (i-2)*(nplaces+ntrans)+ 1 : i*(nplaces+ntrans)) = bloque;
                beq2 = [beq2; zeros(nplaces,1)];
                % add collision avoidance constraints
                Aineq2((i-1)*(nplaces)+1: i*nplaces, (i-2)*(nplaces+ntrans)+1 : i*(nplaces+ntrans)) = bloque2;
                bineq2 = [bineq2; ones(nplaces,1)];
            end
        end
        % put the final marking
        Aeq2(num_intermediate*(nplaces)+1: (num_intermediate+1)*nplaces, ...
            (num_intermediate-1)*(nplaces+ntrans)+ 1 : (num_intermediate-1)*(nplaces+ntrans)+ nplaces) = eye(nplaces);
        beq2 = [beq2;mf];

        tic;
        [sol2,fval2,exitflag2] = intlinprog(f2,[],Aineq2,bineq2,Aeq2,beq2,zeros(1,size(Aeq2,2)),[],[],opt);
        time_LP2 = toc;
        if (exitflag2 ~= 1)
            fprintf(1,'Error solving second LP for %d intermediate markings! Increment the number of intermediate markings.\n',num_intermediate);
            num_intermediate = num_intermediate + 1;
        else
            fprintf('Second LP with intermediate markings problem solved!\n');
        end
    end
    
    if flag_ILP
        tic;
        [sol2i,fval2i,exitflag2] = intlinprog(f2,1:size(Aineq2,2),Aineq2,bineq2,Aeq2,beq2,zeros(1,size(Aeq2,2)),[],[],opt);
        time_IP2 = toc;
        if (exitflag2 ~= 1)
            fprintf(1,'Error solving second ILP for %d intermediate markings! \n',num_intermediate);
            return;
        else
            fprintf('\nSecond ILP with intermediate markings problem solved!\n');
        end
    end

    if flag_ILP 
        % save the data in simulation structure
        simulation.interMILP.runtime = time_IP2;
        simulation.interMILP.cost = 0;
        simulation.interMILP.no_inter_markings = num_intermediate;
    end

    % save the data in simulation structure
    simulation.interLP.runtime = time_LP2;
    simulation.interLP.cost = 0;
    simulation.interLP.no_inter_markings = num_intermediate;

    for i = 1 : num_intermediate
        if flag_ILP
            simulation.interMILP.cost = simulation.interMILP.cost + sum(sol2i((i-1)*(nplaces+ntrans) + 1 + nplaces:i*(nplaces+ntrans)));
        end
        simulation.interLP.cost = simulation.interLP.cost + sum(sol2((i-1)*(nplaces+ntrans) + 1 + nplaces:i*(nplaces+ntrans)));
    end
else % num_intermediare = 1 so no necessary to solve this
    fprintf("\nNumber of intermediary markings is equal with 1. The second problem will not be solved.\n");
    if flag_ILP
        % save the data in simulation structure
        simulation.interMILP.runtime = 0;
        simulation.interMILP.cost = 0;
        simulation.interMILP.no_inter_markings = num_intermediate;
    end

    % save the data in simulation structure
    simulation.interLP.runtime = 0;
    simulation.interLP.cost = 0;
    simulation.interLP.no_inter_markings = num_intermediate;

    if flag_ILP
        simulation.interMILP.cost = simulation.completeLP.cost;
        time_IP2 = 0;
        fval2i=0;
    end
    simulation.interLP.cost = simulation.completeLP.cost;
    time_LP2 = 0;
    fval2=0;
end

if plot_animation
    fprintf(1,'\nSolution LP2 found in %f [secs].',time_LP2);
    fprintf(1,'\nOptimal value LP2 = %s\n',num2str(fval2));

    if flag_ILP
        fprintf(1,'\nSolution IP2 found in %f [secs].',time_IP2);
        fprintf(1,'\nOptimal value IP2 = %s\n',num2str(fval2i));
    end

    for i = 1 : num_intermediate
        if (i == 1)
            current_marking = m0;
        else
            current_marking = round(sol2((nplaces+ntrans)*(i-2)+1: (nplaces+ntrans)*(i-2)+nplaces));
            plot_environment(T.rem_cell(setdiff(find(current_marking),find(mf))),...
                T.rem_cell(find(mf)),T.obstacles,env_limit,T.Vert); %partition
            title(sprintf('Iteration %d',i));
        end
        if num_intermediate == 1
            current_sigma = round(sol1(1:ntrans));
        else
            current_sigma = round(sol2((nplaces+ntrans)*(i-1) + nplaces + 1 : (nplaces+ntrans)*(i-1) + nplaces + ntrans));
        end
        [feasible_sigma, Rob_places, ~, ~] = sigma2trajectories(Pre,Post,current_marking,current_sigma,find(current_marking));
        if ~feasible_sigma
            error('something wrong happened!');
        end
        for j = 1 : length(Rob_places)
            traj = Rob_places{j};
            for k = 1 : length(traj)-1
                plot([T.centr{T.rem_cell(traj(k))}(1) T.centr{T.rem_cell(traj(k+1))}(1)],...
                    [T.centr{T.rem_cell(traj(k))}(2) T.centr{T.rem_cell(traj(k+1))}(2)],'-','LineWidth',1,'Color','magenta');
            end
            if (length(traj) > 1)
                init = T.centr{T.rem_cell(traj(1))};
                fill([init(1)-1.5 init(1)+1.5 init(1)+1.5 init(1)-1.5],[init(2)-1.5 init(2)-1.5 init(2)+1.5 init(2)+1.5],'magenta','EdgeColor','none');
                init = T.centr{T.rem_cell(traj(end))};
                fill([init(1)-1.5 init(1)+1.5 init(1)+1.5 init(1)-1.5],[init(2)-1.5 init(2)-1.5 init(2)+1.5 init(2)+1.5],'magenta','EdgeColor','none');
            end
            %final destinations already reached
            reached = intersect(find(current_marking),find(mf));
            for k=1:length(reached)%final points reached
                fill(T.Vert{T.rem_cell(reached(k))}(1,:),T.Vert{T.rem_cell(reached(k))}(2,:),'green','LineWidth',1);%,'LineStyle','--');,'FaceAlpha',0.5
            end
        end
    end
end