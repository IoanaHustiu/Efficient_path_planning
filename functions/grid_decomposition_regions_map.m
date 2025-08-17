function [C,adj]=grid_decomposition_regions_map(env_bounds,obs_size)

%first begin creating points (vertices) X
% obs_size = 10;
x = (env_bounds(2)-env_bounds(1)) / obs_size;
y = (env_bounds(4)-env_bounds(3)) / obs_size;

C=cell(1,x*y);

for j = 0 : y-1
    for i = 0 : x-1
        C{j*x+i+1} = obs_size*[i i+1 i+1 i; j j j+1 j+1];
        % C{i*x+j+1} = obs_size*[j j+1 j+1 j; i i i+1 i+1];
    end
end

nr_reg = length(C);
adj=sparse(eye(nr_reg)); %self transitions in every cell

for i=1:nr_reg
    for j=i+1:nr_reg %go only through higher indices (adjacency matrix is symmetric)
        if (norm(mean(C{i},2)-mean(C{j},2)) <= 10^10*eps+obs_size) %two common vertices means adjacency
            adj(i,j)=1;
            adj(j,i)=1; %symmetric
        end
    end
end
