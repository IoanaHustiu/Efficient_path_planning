function simulation = solve_LPs_boolean(At,bt,Post,Pre,V,m0,T,env_limit,plot_animation)

% for each simulation the runtime, the cost, the number of iteration/intermediary markings and the cell capacity will be saved
simulation.completeMILP = [];
simulation.completeLP = [];
simulation.interMILP = [];
simulation.interLP = [];

C = Post - Pre;
nplaces = size(Post,1);
ntrans = size(Post,2);
np = size(At,2);
nc = size(At,1);
for i = 1 : length(T.Vert)
    T.centr{i} = mean(T.Vert{i},2)';
end
if (plot_animation)
    plot_environment(T.R0,T.props,T.obstacles,env_limit,T.Vert); %partition
    title(sprintf('Iteration %d',0));
end

fprintf(1,'Solve first problem\n');

%%%%first problem
% f1 = [ones(1,ntrans) sum(m0)*100];
% %state equation
% Aeq1 = [C zeros(nplaces,1)];
% beq1 = mf-m0;
% Aineq1 = [Post -ones(nplaces,1)];
% bineq1 = -m0;

f1 = [zeros(1,nplaces) ones(1,ntrans) zeros(1,np) 100];

Aeq1 = [eye(nplaces) -C zeros(nplaces,np) zeros(nplaces,1)];
beq1 = m0;
Aineq1 = [-V zeros(np,ntrans) eye(np) zeros(np,1);...
    V zeros(np,ntrans) -100*eye(np) zeros(np,1);...
    zeros(nc,nplaces) zeros(nc,ntrans) At zeros(nc,1);...
    zeros(nplaces) Post zeros(nplaces,np) -ones(nplaces,1)];
bineq1 = [zeros(np,1);zeros(np,1);bt;-m0];

opt = optimoptions('intlinprog','Display','none');

tic;
[sol1,fval1,exitflag1] = intlinprog(f1,[],Aineq1,bineq1,Aeq1,beq1,zeros(1,size(Aeq1,2)),[inf(nplaces+ntrans,1);ones(np,1);inf],[],opt);
time_LP1 = toc;
if (exitflag1 ~= 1)
    error('Error solving first problem!\n');
else
    fprintf('First LP problem solved!\n');
end

tic;
[sol1i,fval1i,exitflag1i] = intlinprog(f1,1:size(Aeq1,2),Aineq1,bineq1,Aeq1,beq1,zeros(1,size(Aeq1,2)),[inf(nplaces+ntrans,1);ones(np,1);inf],[],opt);
time_IP1 = toc;

if (exitflag1i ~= 1)
    fprintf('Error solving first IP problem!\n');
    return;
else
    fprintf('First IP problem solved!\n');    
end

% save the data in simulation structure
simulation.completeMILP.runtime = time_IP1;
simulation.completeMILP.cost = sum(sol1i(nplaces+1:nplaces+ntrans));
simulation.completeMILP.sol = sol1i(nplaces+1:nplaces+ntrans);
simulation.completeMILP.cellCapacity = sol1i(end);
simulation.completeLP.runtime = time_LP1;
simulation.completeLP.cost = sum(sol1(nplaces+1:nplaces+ntrans));
simulation.completeLP.sol = sol1(nplaces+1:nplaces+ntrans);
simulation.completeLP.cellCapacity = sol1(end);
simulation.LPfixedS.runtime = [];
simulation.LPfixedS.cost = [];
simulation.LPfixedS.cellCapacity = [];
simulation.LPfixedS.iterations = 0;

%check if the solution of the first LP is integer
% if norm(ceil(sol1(end)-eps*1000)- sol1(end)) > 1000*eps %solution is not integer
%     %fix s
%     fprintf(1,'Solve first problem with ceil(s)\n');
%     Aeq1 = [Aeq1; [zeros(1,ntrans) 1]];
%     beq1 = [beq1; ceil(sol1(end)-eps*1000)];
%     tic;
%     [sol1e,~,~] = intlinprog(f1,[],Aineq1,bineq1,Aeq1,beq1,zeros(1,size(Aeq1,2)),[],[],opt);
%     time_LP1plus = toc;
%     simulation.completeLPfixingB.runtime = time_LP1plus;
%     simulation.completeLPfixingB.cost = sum(sol1e(1:ntrans));
%     simulation.completeLPfixingB.cellCapacity = ceil(sol1(end)-eps*1000);
% end

% x1 = sol1(nplaces+ntrans+1:nplaces+ntrans+np);
% sol1(end)
if norm(ceil(sol1(end)-eps*1000)- sol1(end)) < 1000*eps %solution is integer
    sigmaStar = round(sol1(nplaces+1:nplaces+ntrans));
    iterationsRounding = 0;
    time_LP2 = 0;
else
    %fix s and solve the LP with s constant
    fprintf(1,'Firing vector is not integer. Solve the second problem with ceil(s)\n');
    s = ceil(sol1(end));
    f1 = f1(1:end-1);
    Aineq1 = Aineq1(:,1:end-1);
    bineq1(end-nplaces+1:end) = s*ones(nplaces,1)-m0;
    Aeq1 = Aeq1(:,1:end-1);
    tic;
    [sol1e,~,exitflag1e] = intlinprog(f1,[],Aineq1,bineq1,Aeq1,beq1,zeros(1,size(Aeq1,2)),[inf(nplaces+ntrans,1);ones(np,1)],[],opt);
    time_LP1plus = toc;

    if exitflag1e~=1
        fprintf('Error solving second LP problem!\n');
        return;
    else
        fprintf('Second LP problem solved!\n');
    end

    %if sigma is integer
    if norm(ceil(sol1(nplaces+1:nplaces+ntrans)-eps*1000) - sol1(nplaces+1:nplaces+ntrans)) < 1000*eps %solution is integer

        simulation.LPfixedS.runtime = time_LP1plus;
        simulation.LPfixedS.cost = sum(sol1e(nplaces+1:nplaces+ntrans));
        simulation.LPfixedS.cellCapacity = s;
        simulation.LPfixedS.iterations = 1;

        sigmaStar = round(sol1(nplaces+1:nplaces+ntrans));
        iterationsRounding = 0;
        time_LP2 = time_LP1plus;

    else %solution is not integer (ceil one non-integer element and save all integer elements)
        sigmaStar = zeros(ntrans,1);
        integerPos = [];
        nonIntegerPos = [];
        for i=1:ntrans
            if norm(ceil(sol1(nplaces+i)-eps*1000) - sol1(nplaces+i)) < 1000*eps %element is integer
                integerPos = [integerPos i];
                sigmaStar(i) = ceil(sol1(nplaces+i));
            else
                nonIntegerPos = [nonIntegerPos i];
                sigmaStar(i) = 0;
            end
        end
        %ceil one random element
        randomIndex = randi([1,length(nonIntegerPos)]);
        sigmaStar(nonIntegerPos(randomIndex)) = ceil(sol1(nplaces+nonIntegerPos(randomIndex)));
        nonIntegerPos(randomIndex) = [];
        integerPos = sort([integerPos randomIndex]);

        flagNonIntegerSolution = 1;
        iterationsRounding = 1;
        time_LP2 = time_LP1plus;

        while flagNonIntegerSolution
            %reformulate LP2 problem
            bl1 = ones(ntrans,1);
            bl1(nonIntegerPos) = 0;
            Aeq1 = [Aeq1;zeros(ntrans,nplaces) diag(bl1) zeros(ntrans,np)];
            beq1 = [beq1;sigmaStar];

            tic;
            [sol1f,~,exitflag1f] = intlinprog(f1,[],Aineq1,bineq1,Aeq1,beq1,zeros(1,size(Aeq1,2)),[inf(nplaces+ntrans,1);ones(np,1)],[],opt);
            time_LP2 = time_LP2+toc;

            if exitflag1f~=1
                fprintf('Error solving second LP problem! - in while\n');
                return;
            else
                fprintf('Second LP problem solved!\n');
                iterationsRounding = iterationsRounding+1;
                fprintf('Iteration %i\n',iterationsRounding);
            end

            i=1;
            while i<=length(nonIntegerPos)
                if norm(ceil(sol1f(nplaces+nonIntegerPos(i))-1000*eps) - sol1f(nplaces+nonIntegerPos(i))) < 1000*eps %element is integer
                    integerPos = [integerPos nonIntegerPos(i)];
                    sigmaStar(nonIntegerPos(i)) = ceil(sol1f(nplaces+nonIntegerPos(i)));
                    nonIntegerPos(i) = [];
                else
                    i = i+1;
                end
            end

            if ~isempty(nonIntegerPos)
                %ceil one random element
                randomIndex = randi([1,length(nonIntegerPos)]);
                sigmaStar(nonIntegerPos(randomIndex)) = ceil(sol1f(nplaces+nonIntegerPos(randomIndex)));
                nonIntegerPos(randomIndex) = [];
                integerPos = sort([integerPos randomIndex]);

            else
                flagNonIntegerSolution = 0;
                simulation.LPfixedS.runtime = time_LP2;
                simulation.LPfixedS.cost = sum(sol1f(nplaces+1:nplaces+ntrans));
                simulation.LPfixedS.cellCapacity = s;
                simulation.LPfixedS.iterations = iterationsRounding;
            end
        end
    end
end

simulation.LPfixedS.sol = sigmaStar;
mf = m0 + C*sigmaStar;


if plot_animation
    fprintf(1,'\n=============================================');
    fprintf(1,'\n           Solving using intlinprog');
    fprintf(1,'\n=============================================');
    fprintf(1,'\nSolution LP1 found in %f [secs].',time_LP1);
    fprintf(1,'\nOptimal solution LP1 = %s',num2str(fval1));
    fprintf(1,'\nInfinite norm LP1 of Post * sigma=%s',num2str(sol1(end)));

    fprintf(1,'\n\nSolution IP1 found in %f [secs].',time_IP1);
    fprintf(1,'\nOptimal solution IP1 = %s',num2str(fval1i));
    fprintf(1,'\nInfinite norm IP1 of Post * sigma=%s',num2str(sol1i(end)));

    fprintf(1,'\n\nSolution LP2 found in %f [secs].',time_LP2);
%     fprintf(1,'\nOptimal solution LP2 = %s',num2str(sum(sol1f(nplaces+1:nplaces+ntrans))));
%     fprintf(1,'\nInfinite norm LP2 of Post * sigma=%s',num2str(sol1i(end)));
    fprintf(1,'\nNumber of iterations of LP2 for obtaining integer solution: %i',iterationsRounding)
    
    plot_environment(setdiff(find(m0),find(mf)),T.props,T.obstacles,env_limit,T.Vert); %plot partition
    title(sprintf('Final marking'));

    [feasible_sigma, Rob_places, ~, ~] = sigma2trajectories(Pre,Post,m0,sigmaStar,find(m0));
    if ~feasible_sigma
        error('something wrong happened!');
    end
    for j = 1 : length(Rob_places)
        traj = Rob_places{j};
        for k = 1 : length(traj)-1
            if k==length(traj)-1
                fill(T.Vert{traj(k+1)}(1,:),T.Vert{traj(k+1)}(2,:),'green','LineWidth',1);
            end
            plot([T.centr{traj(k)}(1) T.centr{traj(k+1)}(1)],[T.centr{traj(k)}(2) T.centr{traj(k+1)}(2)],'-','LineWidth',1,'Color','magenta');
        end
        if (length(traj) > 1)
            init = T.centr{traj(1)};
            fill([init(1)-1.5 init(1)+1.5 init(1)+1.5 init(1)-1.5],[init(2)-1.5 init(2)-1.5 init(2)+1.5 init(2)+1.5],'magenta','EdgeColor','none');
            init = T.centr{traj(end)};
            fill([init(1)-1.5 init(1)+1.5 init(1)+1.5 init(1)-1.5],[init(2)-1.5 init(2)-1.5 init(2)+1.5 init(2)+1.5],'magenta','EdgeColor','none');
        end
        %final destinations already reached
        reached = intersect(find(m0),find(mf));
        for k=1:length(reached)%final points reached
            fill(T.Vert{reached(k)}(1,:),T.Vert{reached(k)}(2,:),'green','LineWidth',1);%,'LineStyle','--');,'FaceAlpha',0.5
        end
    end
end