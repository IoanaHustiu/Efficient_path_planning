function [LP1,ILP1] = solve_mILPr(Post,Pre,mf,m0,flag_ILP)
LP1.exitflag = 0;
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

    time_limit = 1200;
    opt = optimoptions('intlinprog','Display','none','MaxTime', time_limit);
    upper = [inf(ntrans,1);nR];

    tic;
    [sol1,fval1,exitflag1] = intlinprog(f1,[],Aineq1,bineq1,Aeq1,beq1,zeros(1,size(Aeq1,2)),upper,[],opt);
    time_LP1 = toc;

    if time_LP1 > time_limit
        fprintf('Solving the LP problem takes more than %i seconds!\n',time_limit);
        LP1.exitflag = 3600;
        LP1.runtime = time_LP1;
    else
        if (exitflag1 ~= 1)
            fprintf('Error solving first problem (not for time limit reasons)!\n');
            LP1.exitflag = exitflag1;
            LP1.runtime = time_LP1;
        else
            fprintf('First LP problem solved! (m, sigma, s - real variables)\n');
            LP1.sol = sol1;
            LP1.cost = fval1;
            LP1.exitflag = exitflag1;
            LP1.runtime = time_LP1;
        end
    end

    if flag_ILP
        tic;
        [sol1i,fval1i,exitflag1i] = intlinprog(f1,1:size(Aineq1,2),Aineq1,bineq1,Aeq1,beq1,zeros(1,size(Aeq1,2)),[],[],opt);
        time_ILP1 = toc;

        if time_ILP1 > time_limit
            fprintf('Solving the ILP problem takes more than %i seconds!\n',time_limit);
            ILP1.exitflag = 3600;
            ILP1.runtime = time_ILP1;
        else
            if (exitflag1 ~= 1)
                fprintf('Error solving first problem!\n');
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
    LP1.exitflag = 10;
    ILP1.exitflag = 10;
    return;
end

end