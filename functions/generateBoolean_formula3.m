function [At,bt] = generateBoolean_formula3(nc,np,nr)
    % Input:
    %   nc - number of clauses (rows)
    %   np - number of variables (columns)
    %   nr - number of robots
    % Output:
    %   At - nc x np constraint matrix
    %   bt - nc x 1 column vector
flag_x = 1;

At = [];

% generate the solution of the problem
while flag_x
    x = randi([0 1],nc,1);
    if sum(x)>0 && sum(x)<=nr
        flag_x = 0;
    end
end

for i=1:nc
    if x(i)
        flag_row = 1;
        while flag_row
            row = randi([-1,0],1,np);
            if any(row)
                At = [At;row];
                flag_row = 0;
            end
        end
    end
end

bt = (sum(At'==1)-1)';

end