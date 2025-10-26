function [MILP,ILP2] = solve_LP_IM_boolean(num_intermediate,C,Post,At,bt,V,m0,flag_ILP)
MILP.exitflag = 0;
ILP2.exitflag = 0;
nplaces = size(C,1);
ntrans = size(C,2);
np = size(At,2);
nc = size(At,1);
nR = sum(m0);
M = 500*nR;

try
    Aeq2 = zeros((num_intermediate)*nplaces,num_intermediate*(nplaces+ntrans)+np);
    Aineq2 = zeros(2*np+num_intermediate*nplaces+nc,num_intermediate*(nplaces+ntrans)+np);
    bineq2 = [];

    f2 = [zeros(1,nplaces) ones(1,ntrans)];
    Aeq2(1:nplaces,1:nplaces+ntrans) = [eye(nplaces) -C];
    beq2= m0;
    bloque_M1=[-eye(nplaces) zeros(nplaces,ntrans) eye(nplaces) -C];

    for i = 2 : num_intermediate %(num_intermediate - i + 2)*10*sum(m0)
        f2 = [f2 zeros(1,nplaces) ones(1,ntrans)];
        % add state equation
        Aeq2((i-1)*(nplaces)+1: i*nplaces, (i-2)*(nplaces+ntrans)+ 1 : i*(nplaces+ntrans)) = bloque_M1;
        beq2 = [beq2; zeros(nplaces,1)];
    end

    f2 = [f2 zeros(1,np)];

    %put capacity constraints
    Aineq2(1:nplaces,1:nplaces+ntrans) = [zeros(nplaces,nplaces) Post];
    bineq2 = [bineq2;-m0+ones(nplaces,1)];

    bloque_X=[eye(nplaces) zeros(nplaces,ntrans+nplaces) Post];
    for i = 2 : num_intermediate
        Aineq2((i-1)*nplaces+1:i*nplaces,...
            (i-2)*(nplaces+ntrans)+1 : i*(nplaces+ntrans)) = bloque_X;
        bineq2 = [bineq2;ones(nplaces,1)];
    end

    %add constraints regarding the final marking and x*
    Aineq2(num_intermediate*nplaces+1:num_intermediate*nplaces+np,(num_intermediate-1)*(nplaces+ntrans)+1:(num_intermediate-1)*(nplaces+ntrans)+nplaces) = -V; %x* <= V^T*x
    Aineq2(num_intermediate*nplaces+1:num_intermediate*nplaces+np,num_intermediate*(nplaces+ntrans)+1:num_intermediate*(nplaces+ntrans)+np) = eye(np);
    bineq2(num_intermediate*nplaces+1:num_intermediate*nplaces+np) = zeros(np,1);

    Aineq2(num_intermediate*nplaces+np+1:num_intermediate*nplaces+np+np,(num_intermediate-1)*(nplaces+ntrans)+1:(num_intermediate-1)*(nplaces+ntrans)+nplaces) = V; %V^T*x <= M*x
    Aineq2(num_intermediate*nplaces+np+1:num_intermediate*nplaces+np+np,num_intermediate*(nplaces+ntrans)+1:num_intermediate*(nplaces+ntrans)+np) = -M*eye(np);
    bineq2(num_intermediate*nplaces+np+1:num_intermediate*nplaces+np+np) =  zeros(np,1);

    % add constraints for Boolean formula
    Aineq2(num_intermediate*nplaces+np+np+1:num_intermediate*nplaces+np+np+nc,num_intermediate*(nplaces+ntrans)+1:num_intermediate*(nplaces+ntrans)+np) = At;
    bineq2(num_intermediate*nplaces+np+np+1:num_intermediate*nplaces+np+np+nc) = bt;


    Aeq2 = sparse(Aeq2);
    Aineq2 = sparse(Aineq2);
    beq2 = sparse(beq2);
    bineq2 = sparse(bineq2);

    time_limit = 1200;
    opt = optimoptions('intlinprog','Display','none','MaxTime',time_limit);

    intcon = num_intermediate*(nplaces+ntrans)+1:num_intermediate*(nplaces+ntrans)+np;
    upper = [];
    for i = 1 : num_intermediate
        upper = [upper inf*ones(1,nplaces+ntrans)];
    end
    upper = [upper ones(1,np)];

    fprintf(1,'MILP with %d equality constraints and %d inequality constraints.\n',size(Aeq2,1),size(Aineq2,1));
    tic;
    [sol2,fval2,exitflag2] = intlinprog(f2,intcon,Aineq2,bineq2,Aeq2,beq2,zeros(1,size(Aeq2,2)),upper,[],opt);
    time_LP2 = toc;

    if time_LP2>=time_limit
        fprintf('Solving the second (MILP) problem takes more than %i seconds!\n',time_limit);
        MILP.exitflag = 3600;
        MILP.runtime = time_LP2;
    else
        if (exitflag2 ~= 1)
            fprintf(1,'Error solving second problem (MILP) with %d intermediate markings! Increment the number of intermediate markings.\n',num_intermediate);
            MILP.exitflag = exitflag2;
            MILP.runtime = time_LP2;
            % return;
        else
            fprintf('MILP with intermediate markings problem solved!\n');
            MILP.sol = sol2;
            MILP.cost = fval2;
            MILP.exitflag = exitflag2;
            MILP.runtime = time_LP2;
        end
    end

    if flag_ILP
        tic;
        [sol2i,fval2i,exitflag2i] = intlinprog(f2,1:size(Aineq2,2),Aineq2,bineq2,Aeq2,beq2,zeros(1,size(Aeq2,2)),upper,[],opt);
        time_IP2 = toc;

        if time_LP2 >= time_limit
            fprintf('Solving the second (ILP) problem takes more than %i seconds!\n',time_limit);
            ILP2.exitflag = 3600;
            ILP2.runtime = time_IP2;
        else
            if (exitflag2i ~= 1)
                fprintf(1,'Error solving second ILP for %d intermediate markings! \n',num_intermediate);
                ILP2.exitflag = exitflag2i;
                ILP2.runtime = time_IP2;
            else
                fprintf('Second ILP with intermediate markings problem solved!\n');
                ILP2.sol = sol2i;
                ILP2.cost = fval2i;
                ILP2.exitflag = exitflag2i;
                ILP2.runtime = time_IP2;
            end
        end
    end
catch
    warning('Out of memory error.');
    ILP2.exitflag = 10;
    MILP.exitflag = 10;
    return;
end

end