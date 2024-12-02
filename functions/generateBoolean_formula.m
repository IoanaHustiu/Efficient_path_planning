function [At,bt] = generateBoolean_formula(nc,np)
    % Input:
    %   nc - number of clauses (rows)
    %   np - number of variables (columns)
    % Output:
    %   At - nc x np matrix satisfying all constraints and ensuring satisfiability

    % Initialize the matrix
    At = zeros(nc,np);

    % Step 1: Ensure all rows have at least one nonzero entry
    for i=1:nc
        % Randomly decide the number of literals in this clause
        num_literals = randi([1,np]); % At least one literal
        literal_indices = randperm(np,num_literals);
        
        % Randomly assign -1 or +1 to selected literals
        At(i,literal_indices) = randsample([-1,1],num_literals,true);
    end

    % Step 2: Ensure no row has all entries equal to 1
    for i=1:nc
        if all(At(i,:) == 1)
            % Flip one literal to -1
            pos_indices = find(At(i,:) == 1);
            At(i,pos_indices(1)) = -1;
        end
    end

    % Step 3: Ensure no column has all entries equal to 0
    for j=1:np
        if all(At(:,j) == 0)
            % Randomly assign -1 or +1 to a row in this column
            row_to_update = randi(nc);
            At(row_to_update,j) = randsample([-1, 1],1);
        end
    end

    % Step 4: Ensure at least one column has no +1
    % Identify all columns with +1
    cols_with_positive = any(At == 1,1);
    % Select a column to remove +1 values from (if all columns have +1)
    if all(cols_with_positive)
        col_to_fix = randi(np);
        At(At(:,col_to_fix) == 1,col_to_fix) = -1; % Flip all +1 to -1
    end

    % Step 5: Ensure satisfiability
    % Create a random satisfying assignment
    satisfying_assignment = randi([0, 1],1,np); % 0 for false, 1 for true

    % Validate the matrix against the satisfying assignment
    for i=1:nc
        % Check if the clause is satisfied under the current assignment
        clause = At(i,:);
        clause_eval = (clause == -1 & satisfying_assignment == 1) | ...
                      (clause == 1 & satisfying_assignment == 0);

        % If clause is unsatisfied, adjust it
        if ~any(clause_eval)
            % Randomly select a variable and ensure it satisfies the clause
            var_to_fix = randi(np);
            if satisfying_assignment(var_to_fix) == 1
                At(i,var_to_fix) = -1; % Add negative literal
            else
                At(i,var_to_fix) = 1;  % Add positive literal
            end
        end
    end
    bt = (sum(At'==1)-1)';
end