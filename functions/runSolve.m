function runSolve()

%% linprog
x = optimvar('x');
y = optimvar('y');
prob = optimproblem;
prob.Objective = -x - y/3;
prob.Constraints.cons1 = x + y <= 2;
prob.Constraints.cons2 = x + y/4 <= 1;

opt = optimoptions('linprog','Algorithm','dual-simplex','Display','none');
sol = solve(prob,'Options',opt);

% problem = prob2struct(prob);
% problem.options.Display = 'none';
% x = linprog(problem.f,problem.Aineq,problem.bineq,problem.Aeq,problem.beq,problem.lb,problem.ub,problem.options);

%% intlinprog
x = optimvar('x','Type','integer');
y = optimvar('y','Type','integer');
prob = optimproblem;
prob.Objective = -x - y/3;
prob.Constraints.cons1 = x + y <= 2;
prob.Constraints.cons2 = x + y/4 <= 1;

opt = optimoptions('intlinprog','Display','none');
intcon = [1 2];
sol = solve(prob,'Options',opt);

% problem = prob2struct(prob);
% problem.options.Display = 'none';
% x = intlinprog(problem.f,problem.intcon,problem.Aineq,problem.bineq,problem.Aeq,problem.beq,problem.lb,problem.ub,problem.options);

end