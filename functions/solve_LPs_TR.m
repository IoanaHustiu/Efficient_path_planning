function [simulation,flag] = solve_LPs_TR(Post,Pre,m0,V,T,env_limit,plot_animation,flag_ILP)

% for each simulation the runtime, the cost, the number of iteration/intermediary markings and the cell capacity will be saved
simulation.completeMILP = [];
simulation.completeLP = [];
flag = 1;
C = Post - Pre;
nplaces = size(Post,1);
ntrans = size(Post,2);
nR = sum(m0);
np = size(V,1);
num_intermediary_markings = ceil(np/nR);

for i = 1 : length(T.Vert)
    T.centr{i} = mean(T.Vert{i},2)';
end
if (plot_animation)
    plot_environment(T.rem_cell(find(m0)),T.props,T.obstacles,env_limit,T.Vert); %partition
    title(sprintf('Iteration %d',1));
end

%% first LP problem

f1 = [];
Aeq1 = zeros(num_intermediary_markings*nplaces,num_intermediary_markings*(nplaces+ntrans)+1);
Aeq1(1:nplaces,1:nplaces+ntrans) = [eye(nplaces) -C];
beq1 = m0;

Aineq1 = zeros(num_intermediary_markings*nplaces+nplaces,num_intermediary_markings*(nplaces+ntrans)+1);
Aineq1(1:nplaces,nplaces+1:nplaces+ntrans) = Post;
Aineq1(1:nplaces,end) = -ones(nplaces,1);
bineq1 = zeros(num_intermediary_markings*nplaces+nplaces,1);
bineq1(1:nplaces) = -m0;
bineq1(num_intermediary_markings*nplaces+1:end) = -V'*ones(np,1);

bloque1 = [-eye(nplaces) zeros(nplaces,ntrans) eye(nplaces) -C];
bloque2 = [eye(nplaces) zeros(nplaces,ntrans) zeros(nplaces) Post];

for i=1:num_intermediary_markings
    %         fprintf("iteration %i\n",i);
    f1 = [f1 zeros(1,nplaces) ones(1,ntrans)];

    if i>1
        Aeq1((i-1)*nplaces+1:i*nplaces,(i-2)*(nplaces+ntrans)+1:i*(nplaces+ntrans)) = bloque1;
        beq1 = [beq1;zeros(nplaces,1)];
    end

    if i>1
        Aineq1((i-1)*nplaces+1:i*nplaces,(i-2)*(nplaces+ntrans)+1:i*(nplaces+ntrans)) = bloque2;
        Aineq1((i-1)*nplaces+1:i*nplaces,end) = -ones(nplaces,1);
    end

    Aineq1(num_intermediary_markings*nplaces+1:end,(i-1)*(nplaces+ntrans)+1:(i-1)*(nplaces+ntrans)+nplaces) = -eye(nplaces);

end

f1 = [f1 100];

opt = optimoptions('intlinprog','Display','none');

tic;
[sol1,fval1,exitflag1,output] = intlinprog(f1,[],Aineq1,bineq1,Aeq1,beq1,zeros(1,size(Aeq1,2)),[],[],opt);
time_LP1 = toc;
if (exitflag1 ~= 1)
    fprintf('Error solving first problem!\n');
else
    fprintf('First LP problem solved!\n');
    simulation.completeLP.runtime = time_LP1;
    simulation.completeLP.cost = sum(f1(1:end-1)*sol1(1:end-1));
    % simulation.completeLP.sol = sol1(nplaces+1:nplaces+ntrans);
    simulation.completeLP.cellCapacity = sol1(end);
    simulation.LPfixedS.runtime = [];
    simulation.LPfixedS.cost = [];
    simulation.LPfixedS.cellCapacity = [];
end

if flag_ILP
    tic;
    [sol1i,fval1i,exitflag1i] = intlinprog(f1,1:size(Aeq1,2),Aineq1,bineq1,Aeq1,beq1,zeros(1,size(Aeq1,2)),[],[],opt);
    time_IP1 = toc;

    if (exitflag1i ~= 1)
        fprintf('Error solving first IP problem!\n');
        return;
    else
        fprintf('First IP problem solved!\n');
        simulation.completeMILP.runtime = time_IP1;
        simulation.completeMILP.cost = sum(f1(1:end-1)*sol1i(1:end-1));
        % simulation.completeMILP.sol = sol1i(nplaces+1:nplaces+ntrans);
        simulation.completeMILP.cellCapacity = sol1i(end);
    end
end


% check if the solution of the first LP is integer
if norm(ceil(sol1(end)-eps*1000)- sol1(end)) > 1000*eps %solution is not integer
    %fix s
    s = ceil(sol1(end)-eps*1000);
    fprintf(1,'Solve second LP problem (first problem with ceil(s))\n');
    f1(end) = [];
    Aeq1(:,end) = [];
    Aineq1(:,end) = [];
    bineq1(1:num_intermediary_markings*nplaces) = bineq1(1:num_intermediary_markings*nplaces)+ceil(sol1(end)-eps*1000);
    tic;
    [sol1e,~,exitflag1e,~] = intlinprog(f1,[],Aineq1,bineq1,Aeq1,beq1,zeros(1,size(Aeq1,2)),[],[],opt);
    time_LP1plus = toc;

    if exitflag1e~=1
        fprintf('Error solving the second LP!\n');
    else
        fprintf('Second LP solved!\n');
        simulation.completeLPfixingB.runtime = time_LP1plus;
        simulation.completeLPfixingB.cost = sum(f1*sol1e);
        simulation.completeLPfixingB.cellCapacity = ceil(sol1(end)-eps*1000);
        sol1 = sol1e;
    end
else
    s = ceil(sol1(end)-eps*1000);
    fprintf('Capacity is integer\n');
end


if plot_animation
    fprintf(1,'\nSolution LP1 found in %f [secs].',time_LP1);
    fprintf(1,'\nOptimal solution LP1 = %s',num2str(fval1));
    fprintf(1,'\nInfinite norm LP1 of Post * sigma=%s',num2str(s));

    if flag_ILP
        fprintf(1,'\n\nSolution IP1 found in %f [secs].',time_IP1);
        fprintf(1,'\nOptimal solution IP1 = %s',num2str(fval1i));
        fprintf(1,'\nInfinite norm IP1 of Post * sigma=%s\n',num2str(sol1i(end)));
    end

    reached = [];
    for i = 1 : num_intermediary_markings
        if (i == 1)
            current_marking = m0;
        else
            current_marking = round(sol1i((nplaces+ntrans)*(i-2)+1:(nplaces+ntrans)*(i-2)+nplaces));
            plot_environment(T.rem_cell(setdiff(find(current_marking),T.props)),...
                T.props,T.obstacles,[0,env_limit,0,env_limit],T.Vert); %partition
            title(sprintf('Iteration %d',i+1));
        end

        if num_intermediary_markings == 1
            current_sigma = round(sol1(nplaces+1:nplaces+ntrans));
        else
            current_sigma = round(sol1((nplaces+ntrans)*(i-1)+nplaces+1:(nplaces+ntrans)*(i-1)+nplaces+ntrans));
        end
        [feasible_sigma,Rob_places,~,~] = sigma2trajectories(Pre,Post,current_marking,current_sigma,find(current_marking));
        if ~feasible_sigma
            error('something wrong happened!');
        end
        %final destinations already reached
        reached = [reached;intersect(find(current_marking),find(sum(V)))];
        for k=1:length(reached)%final points reached
            fill(T.Vert{T.rem_cell(reached(k))}(1,:),T.Vert{T.rem_cell(reached(k))}(2,:),'green','LineWidth',1);%,'LineStyle','--');,'FaceAlpha',0.5
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
            else
                init = T.centr{T.rem_cell(traj(1))};
                fill([init(1)-1.5 init(1)+1.5 init(1)+1.5 init(1)-1.5],[init(2)-1.5 init(2)-1.5 init(2)+1.5 init(2)+1.5],'magenta','EdgeColor','none');
            end
        end
    end
end
end