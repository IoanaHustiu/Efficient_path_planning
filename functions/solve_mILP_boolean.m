function [MILP1,ILP1] = solve_mILP_boolean(Post,Pre,V,At,bt,m0,flag_ILP)
MILP1.exitflag = 0;
ILP1.exitflag = 0;

Post = sparse(Post);
Pre = sparse(Pre);

try
    C = Post - Pre;
    nplaces = size(Post,1);
    ntrans = size(Post,2);
    nc = size(At,1);
    np = size(At,2);
    nR = sum(m0);
    M = 500*nR;

    %construct matrices
    f1 = [zeros(nplaces,1);ones(ntrans,1);zeros(np,1);M];
    Aeq1 = sparse([eye(nplaces) -C zeros(nplaces,np) zeros(nplaces,1)]);
    beq1 = sparse(m0);
    Aineq1 = sparse([-V zeros(np,ntrans) eye(np) zeros(np,1);...
        V zeros(np,ntrans) -M*eye(np) zeros(np,1);...
        zeros(nc,nplaces) zeros(nc,ntrans) At zeros(nc,1);...
        zeros(nplaces) Post zeros(nplaces,np) -ones(nplaces,1)]);
    bineq1 = sparse([zeros(np,1);zeros(np,1);bt;-m0]);

    opt = optimoptions('intlinprog','Display','none','MaxTime', 1800);

    % intcon = nplaces+ntrans+1:nplaces+ntrans+np+1;
    upper = [nR*ones(nplaces,1);inf(ntrans,1);ones(np,1);nR];
    fprintf(1,'\nStart solving the LP formulation...\n');
    tic;
    [sol1,fval1,exitflag1] = intlinprog(f1,[],Aineq1,bineq1,Aeq1,beq1,zeros(1,size(Aeq1,2)),upper,[],opt);
    time_LP1 = toc;
    if (exitflag1 ~= 1)
        if time_LP1 >= 1800
            fprintf('Solving the first (LP) problem lasts more than 1800 seconds!\n');
            MILP1.exitflag = 3600;
        else
            fprintf('Error solving the first problem!\n');
        end
    else
        fprintf('First LP problem solved! (m, sigma,, x, s - real variables)\n');
        MILP1.sol = sol1;
        MILP1.cost = fval1;
        MILP1.exitflag = exitflag1;
        MILP1.runtime = time_LP1;
    end


    %% is still necessary to solve (5) - ILP?
    if flag_ILP
        fprintf(1,'\nStart solving the ILP formulation...\n');
        tic;
        [sol1i,fval1i,exitflag1i] = intlinprog(f1,1:size(Aineq1,2),Aineq1,bineq1,Aeq1,beq1,zeros(1,size(Aeq1,2)),[],[],opt);
        time_ILP1 = toc;
        if (exitflag1 ~= 1)
            if time_ILP1 >= 1800
                fprintf('Solving the ILP problem lasts more than 1800 seconds!\n');
                ILP1.exitflag = 3600;
            else
                fprintf('Error solving the first problem - time limit not reached!\n');
            end
        else
            fprintf('ILP problem solved! (m, sigma, x and s - integer variables)\n');
            ILP1.sol = sol1i;
            ILP1.cost = fval1i;
            ILP1.exitflag = exitflag1i;
            ILP1.runtime = time_ILP1;
        end
    else
        ILP1 = 0;
    end
catch
    warning('Out of memory error.');
    ILP1.exitflag = 10;
    MILP1.exitflag = 10;
    return;
end

end