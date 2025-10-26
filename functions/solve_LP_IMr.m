function [LP2,ILP2] = solve_LP_IMr(num_intermediate,C,Post,m0,mf,flag_ILP)
ILP2.exitflag = 0;
LP2.exitflag = 0;
nplaces = size(C,1);
ntrans = size(C,2);
nR = sum(m0);

try
    Aeq2 = zeros((num_intermediate+1)*nplaces,num_intermediate*(nplaces+ntrans));
    Aineq2 = zeros(num_intermediate*nplaces,num_intermediate*(nplaces+ntrans));

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
    % put the final marking
    Aeq2(num_intermediate*(nplaces)+1: (num_intermediate+1)*nplaces, ...
        (num_intermediate-1)*(nplaces+ntrans)+ 1 : (num_intermediate-1)*(nplaces+ntrans)+ nplaces) = eye(nplaces);
    beq2 = [beq2;mf];

    Aeq2 = sparse(Aeq2);
    beq2 = sparse(beq2);
    Aineq2 = sparse(Aineq2);
    bineq2 = sparse(bineq2);

    time_limit = 1200;
    opt = optimoptions('intlinprog','Display','none','MaxTime',time_limit);
    upper = [];
    for i=1:num_intermediate
        upper = [upper nR*ones(1,nplaces) inf*ones(1,ntrans)];
    end

    tic;
    [sol2,fval2,exitflag2] = intlinprog(f2,[],Aineq2,bineq2,Aeq2,beq2,zeros(1,size(Aeq2,2)),[],[],opt);
    time_LP2 = toc;

    if time_LP2 > time_limit
        fprintf('Solving the second (LP) problem takes more than %i seconds!\n',time_limit);
        LP2.exitflag = 3600;
        LP2.runtime = time_LP2;
    else
        if (exitflag2 ~= 1)
            fprintf(1,'Error solving second LP with %d intermediate markings! Increment the number of intermediate markings.\n',num_intermediate);
            LP2.exitflag = exitflag2;
            LP2.runtime = time_LP2;
        else
            fprintf('Second LP with intermediate markings problem solved!\n');
            LP2.sol = sol2;
            LP2.cost = fval2;
            LP2.exitflag = exitflag2;
            LP2.runtime = time_LP2;
        end
    end

    if flag_ILP
        tic;
        [sol2i,fval2i,exitflag2i] = intlinprog(f2,1:size(Aineq2,2),Aineq2,bineq2,Aeq2,beq2,zeros(1,size(Aeq2,2)),[],[],opt);
        time_ILP2 = toc;

        if time_ILP2 > time_limit
            fprintf('Solving the ILP problem (with intermediary markings) takes more than %i seconds!\n',time_limit);
            ILP2.exitflag = 3600;
            ILP2.runtime = time_ILP2;
        else
            if (exitflag2i ~= 1)
                fprintf(1,'Error solving ILP for %d intermediate markings! Increment the number of intermediate markings.\n',num_intermediate);
                ILP2.exitflag = exitflag2i;
                ILP2.runtime = time_ILP2;
            else
                fprintf('ILP with intermediate markings problem solved!\n');
                ILP2.sol = sol2i;
                ILP2.cost = fval2i;
                ILP2.exitflag = exitflag2i;
                ILP2.runtime = time_ILP2;
            end
        end
    end
catch
    warning('Out of memory error.');
    ILP2.exitflag = 10;
    LP2.exitflag = 10;
    return;
end

end