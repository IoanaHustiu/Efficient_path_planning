function [simulation,flag] = solve_LPs_boolean4_randomX(At,bt,Post,Pre,V,m0,T,env_limit,plot_animation,flag_ILP,flag_SP)

% for each simulation the runtime, the cost, the number of iteration/intermediary markings and the cell capacity will be saved
simulation.completeMILP = [];
simulation.completeLP = [];
flag = 1;
C = Post - Pre;
nplaces = size(Post,1);
ntrans = size(Post,2);
np = size(At,2);
nc = size(At,1);
nR = sum(m0);

for i = 1 : length(T.Vert)
    T.centr{i} = mean(T.Vert{i},2)';
end
% if (plot_animation)
%     plot_environment(T.rem_cell(find(m0)),T.props,T.obstacles,env_limit,T.Vert); %partition
%     title(sprintf('Iteration %d',0));
% end

%% first problem
fprintf(1,'Solve first problem\n');
flagx = 1;
nonIntegerPos = [];
integerPos = [];

f1 = [zeros(nplaces,1);ones(ntrans,1);zeros(np,1);100];
Aeq1 = [eye(nplaces) -C zeros(nplaces,np) zeros(nplaces,1);...
    zeros(np,nplaces) zeros(np,ntrans) zeros(np) zeros(np,1)];
beq1 = [m0;zeros(np,1)];
Aineq1 = [-V zeros(np,ntrans) eye(np) zeros(np,1);...
    V zeros(np,ntrans) -100*eye(np) zeros(np,1);...
    zeros(nc,nplaces) zeros(nc,ntrans) At zeros(nc,1);...
    zeros(nplaces) Post zeros(nplaces,np) -ones(nplaces,1)];
bineq1 = [zeros(np,1);zeros(np,1);bt;-m0];


opt = optimoptions('intlinprog','Display','none');

tic;
[sol1,~,exitflag1] = intlinprog(f1,[],Aineq1,bineq1,Aeq1,beq1,zeros(1,size(Aeq1,2)),[nR*ones(nplaces,1);inf(ntrans,1);ones(np,1);nR],[],opt);
time_LP1 = toc;
if (exitflag1 ~= 1)
    fprintf('Error solving first problem!\n');
    flag = 0;
else
    fprintf('First LP problem solved!\n');
end

if flag_ILP
    tic;
    [sol1i,~,exitflag1i] = intlinprog(f1,1:size(Aeq1,2),Aineq1,bineq1,Aeq1,beq1,zeros(1,size(Aeq1,2)),[nR*ones(nplaces,1);inf(ntrans,1);ones(np,1);nR],[],opt);
    time_IP1 = toc;

    if (exitflag1i ~= 1)
        fprintf('Error solving first IP problem!\n');
    else
        fprintf('First IP problem solved!\n');
    end
end

time_LP1x = 0;
time_LP1xr = 0;
iterationsX = 0;

x = sol1(nplaces+ntrans+1:nplaces+ntrans+np);

while flagx
    iterationsX = iterationsX + 1;
    if norm(ceil(x-eps*1000)-x) > 1000*eps %x is not integer
        for i=1:np %save all integer elements
            if norm((x(i)-eps*1000)-1) < 10000*eps
                integerPos = [integerPos i];
            else
                nonIntegerPos = [nonIntegerPos i];
            end
        end
        nonIntegerPosMax = nonIntegerPos(x(nonIntegerPos)==max(x(nonIntegerPos)));
        if ~isempty(nonIntegerPosMax) %pick one non-integer element (the closest to 1)
            idx = randi(length(nonIntegerPosMax));
            x(nonIntegerPosMax(idx)) = ceil(x(nonIntegerPosMax(idx)));
            nonIntegerPos(find(nonIntegerPos==nonIntegerPosMax(idx))) = [];
            integerPos = [integerPos nonIntegerPosMax(idx)];
            nonIntegerPosMax(idx) = [];

            xStar = zeros(np,1);
            xStar(integerPos) = 1;
            bloque1 = zeros(np,1);
            bloque1(integerPos) = 1;
            bloque1 = diag(bloque1);
            Aeq1(nplaces+1:nplaces+np,:) = [zeros(np,nplaces) zeros(np,ntrans) bloque1 zeros(np,1)];
            beq1(nplaces+1:nplaces+np) = xStar;

            tic;
            [sol1r,~,exitflag1r] = intlinprog(f1,[],Aineq1,bineq1,Aeq1,beq1,zeros(1,size(Aeq1,2)),[nR*ones(nplaces,1);inf(ntrans,1);ones(np,1);nR],[],opt);
            time_LP1xr = time_LP1xr + toc;

            if exitflag1r~=1
                fprintf('Error when rounding x at iteration %i\n',iterationsX);
            else
                if norm((ceil(sol1r(nplaces+ntrans+1:nplaces+ntrans+np)-eps*1000))-sol1r(nplaces+ntrans+1:nplaces+ntrans+np)) < 1000*eps
                    xStar = ceil(sol1r(nplaces+ntrans+1:nplaces+ntrans+np));
                    index = xStar~=0;
                    xStar(index) = 1;
                    s = ceil(sol1r(end)-eps*1000); %fix s
                    flagx = 0;
                else
                    integerPos = [];
                    nonIntegerPos = [];
                end
            end
        end
    else
        xStar = ceil(x-eps*1000);
        s = ceil(sol1(end)-eps*1000); %fix s
        flagx = 0;
    end

end

fprintf(1,'\n=============================================');
fprintf(1,'\n           Reachability case');
fprintf(1,'\n=============================================\n');

clear f1 Aeq1 beq1 Aineq1 beq1


f1 = [zeros(nplaces,1);ones(ntrans,1)];
Aeq1 = [eye(nplaces) -C];
beq1 = m0;
Aineq1 = [zeros(nplaces) Post;...
    -eye(nplaces) zeros(nplaces,ntrans)];
bineq1 = [s-m0;-V'*xStar];

opt = optimoptions('intlinprog','Display','none');

tic;
[sol1x,~,exitflag1x] = intlinprog(f1,[],Aineq1,bineq1,Aeq1,beq1,zeros(1,size(Aeq1,2)),[nR*ones(nplaces,1);inf(ntrans,1)],[],opt);
time_LP2 = toc;
if (exitflag1x ~= 1)
    fprintf('Error solving second problem (x and s are fixed)!\n');
else
    fprintf('Second LP problem solved (x and s are fixed)!\n');
end

% save the data in simulation structure
if flag_ILP
    simulation.completeMILP.runtime = time_IP1;
    simulation.completeMILP.cost = sum(sol1i(nplaces+1:nplaces+ntrans));
    simulation.completeMILP.sol = sol1i(nplaces+1:nplaces+ntrans);
    simulation.completeMILP.cellCapacity = sol1i(end);
    simulation.completeMILP.x = sol1i(nplaces+ntrans+1:nplaces+ntrans+np);
end
simulation.completeLP.runtime = time_LP1;
simulation.completeLP.cost = sum(sol1(nplaces+1:nplaces+ntrans));
simulation.completeLP.sol = sol1(nplaces+1:nplaces+ntrans);
simulation.completeLP.cellCapacity = sol1(end);
simulation.completeLPfixingB.runtime = time_LP1xr;
simulation.completeLPfixingB.runtimeLP2 = time_LP2;
simulation.completeLPfixingB.cost = sum(sol1x(nplaces+1:nplaces+ntrans));
simulation.completeLPfixingB.sol = sol1x(nplaces+1:nplaces+ntrans);
simulation.completeLPfixingB.cellCapacity = s;
simulation.completeLPfixingB.iterations = iterationsX;

fprintf("A number of %i robots will move in order to fulfill the Boolean specification.\n",sum(xStar));

if plot_animation
    fprintf(1,'\n============================================================');
    fprintf(1,'\n    Solution for LP obtained using intlinprog solver ');
    fprintf(1,'\n============================================================');
    fprintf(1,'\nSolution LP1 found in %f [secs].',time_LP1);

    if flag_ILP
        fprintf(1,'\n\nSolution IP1 found in %f [secs].',time_IP1);
        fprintf(1,'\nOptimal solution IP1 = %s',num2str(sum(sol1i(nplaces+1:nplaces+ntrans))));
        fprintf(1,'\nInfinite norm IP1 of Post * sigma = %s',num2str(sol1i(end)));
    end

    fprintf(1,'\n\nSolution LP2 found in %f [secs].',time_LP1xr);
    fprintf(1,'\nOptimal solution LP2 = %s',num2str(sum(sol1x(nplaces+1:nplaces+ntrans))));
    fprintf(1,'\nCell capacity LP2 = %s\n',num2str(s));

    if flag_ILP
        if sum(round(sol1i(nplaces+1:nplaces+ntrans)))>sum(sol1x(nplaces+1:nplaces+ntrans))
            fprintf(2,'\nLP cost is smaller!\n');
        end
    end

    reached = [];
    current_marking = m0;

    if flag_SP
        plot_environment_smartPlant(T.nj,T.rem_cell(setdiff(find(current_marking),T.props)),...
            T.props,T.obstacles,env_limit,T.Vert); %partition
    else
        plot_environment(T.rem_cell(setdiff(find(current_marking),T.props)),...
            T.props,T.obstacles,env_limit,T.Vert); %partition
    end
    title(sprintf('LP solution'));

    current_sigma = round(sol1x(nplaces+1:nplaces+ntrans));

    [feasible_sigma,Rob_places,~,~] = sigma2trajectories(Pre,Post,current_marking,current_sigma,find(current_marking));
    if ~feasible_sigma
        error('something wrong happened!');
    end
    %final destinations already reached
    reached = [reached;intersect(find(current_marking),find(sol1x(1:nplaces)))];
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

if plot_animation && flag_ILP
    reached = [];
    current_marking = m0;

    if flag_SP
        plot_environment_smartPlant(T.nj,T.rem_cell(setdiff(find(current_marking),T.props)),...
            T.props,T.obstacles,env_limit,T.Vert); %partition
    else
        plot_environment(T.rem_cell(setdiff(find(current_marking),T.props)),...
            T.props,T.obstacles,env_limit,T.Vert); %partition
    end
    title(sprintf('ILP solution'));

    current_sigma = round(sol1i(nplaces+1:nplaces+ntrans));

    [feasible_sigma,Rob_places,~,~] = sigma2trajectories(Pre,Post,current_marking,current_sigma,find(current_marking));
    if ~feasible_sigma
        error('something wrong happened!');
    end
    %final destinations already reached
    reached = [reached;intersect(find(current_marking),find(sol1i(1:nplaces)))];
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