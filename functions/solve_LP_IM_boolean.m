function [LPIM,ILPIM] = solve_LP_IM_boolean(num_intermediate,C,Post,At,V,x,m0,flag_ILP)
LPIM.exitflag = 0;
ILPIM.exitflag = 0;
nplaces = size(C,1);
ntrans = size(C,2);
np = size(At,2);
nR = sum(m0);
M = 500*nR;

try
    Aeq2 = zeros((num_intermediate)*nplaces+np,num_intermediate*(nplaces+ntrans+np));
    Aineq2 = zeros(2*num_intermediate*np+num_intermediate*nplaces,num_intermediate*(nplaces+ntrans));

    f2 = [zeros(1,nplaces) ones(1,ntrans) zeros(1,np)];
    Aeq2(1:nplaces,1:nplaces+ntrans) = [eye(nplaces) -C];
    beq2= m0;
    bloque_V1 = [V zeros(np,ntrans) -M*eye(np)]; % Vm - Mx <= 0
    Aineq2(1:np,1:nplaces+ntrans+np) = bloque_V1;
    bineq2 = zeros(np,1);
    bloque_V2 = [-V zeros(np,ntrans) eye(np)]; % -Vm + x <= 0
    Aineq2(np+1:2*np,1:nplaces+ntrans+np) = bloque_V2;
    bineq2 = [bineq2; zeros(np,1)];

    bloque_M1=[-eye(nplaces) zeros(nplaces,ntrans) zeros(nplaces,np) eye(nplaces) -C zeros(nplaces,np)];

    for i = 2 : num_intermediate %(num_intermediate - i + 2)*10*sum(m0)*
        f2 = [f2 zeros(1,nplaces) ones(1,ntrans) zeros(1,np)];
        % add state equation
        Aeq2((i-1)*(nplaces)+1: i*nplaces, (i-2)*(nplaces+ntrans+np)+ 1 : i*(nplaces+ntrans+np)) = bloque_M1;
        beq2 = [beq2; zeros(nplaces,1)];
        % add relations m-x constraints
        Aineq2(2*(i-1)*(np)+1: 2*(i-1)*(np)+np, (i-1)*(nplaces+ntrans+np)+1 : (i)*(nplaces+ntrans+np)) = bloque_V1;
        bineq2 = [bineq2; zeros(np,1)];
        Aineq2(2*(i-1)*(np)+np+1: 2*(i-1)*(np)+2*np, (i-1)*(nplaces+ntrans+np)+1 : (i)*(nplaces+ntrans+np)) = bloque_V2;
        bineq2 = [bineq2; zeros(np,1)];
    end

    %put the Boolean formula
    bloque_X = [zeros(np,nplaces+ntrans) eye(np)];
    constr_x = bloque_X;
    for i = 2 : num_intermediate
        constr_x = [constr_x bloque_X];
    end
    Aeq2(num_intermediate*nplaces+1:num_intermediate*nplaces+np,:) = constr_x;
    beq2 = [beq2;x];
    % Aineq2(num_intermediate*2*np+1:num_intermediate*2*np+np,:) = -constr_x;
    % bineq2 = [bineq2;-x];

    %put capacity constraints
    Aineq2(num_intermediate*2*np+1:num_intermediate*2*np+nplaces,...
        1:nplaces+ntrans) = [zeros(nplaces,nplaces) Post];
    bineq2 = [bineq2;-m0+ones(nplaces,1)];

    bloque_X=[eye(nplaces) zeros(nplaces,ntrans+nplaces+np) Post zeros(nplaces,np)];
    for i = 2 : num_intermediate
        Aineq2((num_intermediate*2)*np+(i-1)*nplaces+1:...
            (num_intermediate*2)*np+(i-1)*nplaces+nplaces,...
            (i-2)*(nplaces+ntrans+np)+1 : i*(nplaces+ntrans+np)) = bloque_X;
        bineq2 = [bineq2;ones(nplaces,1)];
    end

    Aeq2 = sparse(Aeq2);
    Aineq2 = sparse(Aineq2);
    beq2 = sparse(beq2);
    bineq2 = sparse(bineq2);

    opt = optimoptions('intlinprog','Display','none','MaxTime',3600);

    intcon = [];
    upper = [];
    for i = 1 : num_intermediate
        intcon = [intcon (i-1)*(nplaces+ntrans+np)+(nplaces+ntrans)+1:(i-1)*(nplaces+ntrans+np)+(nplaces+ntrans)+np];
        upper = [upper [inf*ones(1,nplaces+ntrans) ones(1,np)]];
    end

    fprintf(1,'MILP with %d variables (%d integer), %d equality constraints and %d inequality constraints.\n',size(Aeq2,2),length(intcon),size(Aeq2,1),size(Aineq2,1));
    tic;
    [sol2,fval2,exitflag2] = intlinprog(f2,intcon,Aineq2,bineq2,Aeq2,beq2,zeros(1,size(Aeq2,2)),upper,[],opt);
    time_LP2 = toc;

    if (exitflag2 ~= 1)
        fprintf(1,'Error solving second MILP for %d intermediate markings! Increment the number of intermediate markings.\n',num_intermediate);
        LPIM.exitflag = 0;
        return;
    else
        if time_LP2>=3600
            fprintf('Solving the first (MILP) problem lasts more than 3600 seconds!\n');
            LPIM.exitflag = 3600;
        else
            fprintf('Second MILP with intermediate markings problem solved!\n');
            LPIM.sol = sol2;
            LPIM.cost = fval2;
            LPIM.exitflag = exitflag2;
            LPIM.runtime = time_LP2;
        end
    end

    if flag_ILP
        tic;
        [sol2i,fval2i,exitflag2i] = intlinprog(f2,1:size(Aineq2,2),Aineq2,bineq2,Aeq2,beq2,zeros(1,size(Aeq2,2)),upper,[],opt);
        time_IP2 = toc;
        if (exitflag2 ~= 1)
            fprintf(1,'Error solving second ILP for %d intermediate markings! \n',num_intermediate);
            ILPIM.exitflag = 0;
            % return;
        else
            if time_LP2>=3600
                fprintf('Solving the first (MILP) problem lasts more than 3600 seconds!\n');
                ILPIM.exitflag = 3600;
            else
                fprintf('Second ILP with intermediate markings problem solved!\n');
                ILPIM.sol = sol2i;
                ILPIM.cost = fval2i;
                ILPIM.exitflag = exitflag2i;
                ILPIM.runtime = time_IP2;
            end
        end
    end
catch
    warning('Out of memory error.');
    ILPIM.exitflag = 10;
    LPIM.exitflag = 10;
    return;
end

end