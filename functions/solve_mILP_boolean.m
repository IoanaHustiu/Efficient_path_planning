function [LP1,ILP1] = solve_mILP_boolean(Post,Pre,V,At,bt,m0,flag_ILP)
LP1.exitflag = 0;
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

    time_limit = 1200;
    opt = optimoptions('intlinprog','Display','none','MaxTime', time_limit);

    % intcon = nplaces+ntrans+1:nplaces+ntrans+np+1;
    upper = [nR*ones(nplaces,1);inf(ntrans,1);ones(np,1);nR];
    fprintf(1,'\nStart solving the LP formulation...\n');
    tic;
    [sol1,fval1,exitflag1] = intlinprog(f1,[],Aineq1,bineq1,Aeq1,beq1,zeros(1,size(Aeq1,2)),upper,[],opt);
    time_LP1 = toc;

    if time_LP1 >= time_limit
        fprintf('Solving the first (LP) problem takes more than %i seconds!\n',time_limit);
        LP1.exitflag = 3600;
        LP1.runtime = time_LP1;
    else
        if (exitflag1 ~= 1)
            fprintf('Error solving the first problem (not for time limit reasons)!\n');
            LP1.exitflag = exitflag1;
            LP1.runtime = time_LP1;
        else
            fprintf('First LP problem solved! (m, sigma, x, s - real variables)\n');
            LP1.sol = sol1;
            LP1.cost = fval1;
            LP1.exitflag = exitflag1;
            LP1.runtime = time_LP1;
        end
    end

    %% solve (5) - ILP
    if flag_ILP
        fprintf(1,'\nStart solving the ILP formulation...\n');
        tic;
        [sol1i,fval1i,exitflag1i] = intlinprog(f1,1:size(Aineq1,2),Aineq1,bineq1,Aeq1,beq1,zeros(1,size(Aeq1,2)),upper,[],opt);
        time_ILP1 = toc;

        if time_ILP1 >= time_limit
            fprintf('Solving the ILP problem takes more than %i seconds!\n',time_limit);
            ILP1.exitflag = 3600;
            ILP1.runtime = time_ILP1;
        else
            if (exitflag1i ~= 1)
                fprintf('Error solving the first problem (not for time limit reasons)!\n');
                ILP1.exitflag = exitflag1i;
                ILP1.runtime = time_ILP1;
            else
                fprintf('ILP problem solved! (m, sigma, x and s - integer variables)\n');
                ILP1.sol = sol1i;
                ILP1.cost = fval1i;
                ILP1.exitflag = exitflag1i;
                ILP1.runtime = time_ILP1;
            end
        end
    end
catch
    warning('Out of memory error.');
    ILP1.exitflag = 10;
    LP1.exitflag = 10;
    return;
end

end