function [MILP1,ILP1] = solve_mILPr(Post,Pre,mf,m0,flag_ILP)
MILP1.exitflag = 0;
ILP1.exitflag = 0;

%initializations
C = Post - Pre;
nplaces = size(Post,1);
ntrans = size(Post,2);
nR = sum(m0);

try
    f1 = [ones(1,ntrans) sum(m0)*100*nR];
    Aeq1 = sparse([C zeros(nplaces,1)]);
    beq1 = sparse(mf-m0);
    Aineq1 = sparse([Post -ones(nplaces,1)]);
    bineq1 = sparse(-m0);

    opt = optimoptions('intlinprog','Display','none','MaxTime', 3600);
    upper = [inf(ntrans,1);nR];

    tic;
    [sol1,fval1,exitflag1] = intlinprog(f1,[],Aineq1,bineq1,Aeq1,beq1,zeros(1,size(Aeq1,2)),upper,[],opt);
    time_LP1 = toc;
    if (exitflag1 ~= 1)
        fprintf('Error solving first problem!\n');
    else
        if time_LP1>3600
            fprintf('Solving the LP problem lasts more than 3600 seconds!\n');
            MILP1.exitflag = 3600;
        else
            fprintf('First LP problem solved! (m, sigma, s - real variables)\n');
            MILP1.sol = sol1;
            MILP1.cost = fval1;
            MILP1.exitflag = exitflag1;
            MILP1.runtime = time_LP1;
        end
    end

    if flag_ILP
        tic;
        [sol1i,fval1i,exitflag1i] = intlinprog(f1,1:size(Aineq1,2),Aineq1,bineq1,Aeq1,beq1,zeros(1,size(Aeq1,2)),[],[],opt);
        time_ILP1 = toc;
        if (exitflag1 ~= 1)
            fprintf('Error solving first problem!\n');
        else
            if time_ILP1>3600
                fprintf('Solving the ILP problem lasts more than 3600 seconds!\n');
                ILP1.exitflag = 3600;
            else
                fprintf('ILP problem solved! (m, sigma, x and s - integer variables)\n');
                ILP1.sol = sol1i;
                ILP1.cost = fval1i;
                ILP1.exitflag = exitflag1i;
                ILP1.runtime = time_ILP1;
            end
        end
    else
        ILP1 = 0;
    end
catch
    warning('Out of memory error.');
    MILP1.exitflag = 10;
    ILP1.exitflag = 10;
    return;
end

end