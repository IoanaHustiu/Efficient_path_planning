function obstacles = random_env_generator(N_o,N_regions)
% Function to create Random Enviroment


obstacles = [];
% generate obstacles in the environment

while (length(obstacles) < N_o)
    obstac = unique(randi([1 N_regions], 1, N_o-length(obstacles)),'stable');
    obstacles = unique([obstacles obstac],'stable');
end

return


adj(obstacles,:) = 0;
adj(:,obstacles) = 0;
% find connected components
G = graph(adj);
[str_bins, binsize] = conncomp(G);
I = binsize(str_bins) == max(binsize);

obstacles = find(I~=1);
if floor(max(binsize)/2) < N_o
    fprintf('Not enough connected cells. New environment is generated.\n');
    obstacles = [];
    N_o = abs(N_o - 10);
else
    flag = 0;
end

% generate initial positions for the robots
    while length(random) < N_o+rob_no
        robots = randi([1 round(0.5*length(rem_cells))], 1, N_o+rob_no-length(random));
        random = unique([random rem_cells(robots)],'stable');  
    end

% generate regions of interest for the robots
    while length(random) < N_o+rob_no+obs_no
        reg = randi([round(0.5*length(rem_cells))+1 length(rem_cells)], 1, N_o+rob_no+obs_no-length(random));
        random = unique([random rem_cells(reg)],'stable');  
    end


for i=1:N_o %obstacles
    cel = random(i);
    pos_x = mod(cel-1,x);
    pos_y = floor((cel-1)/x);
    obstacles{i}=obs_size*[pos_x pos_x+1 pos_x+1 pos_x; pos_y pos_y pos_y+1 pos_y+1];
    o_cells = [o_cells random(i)];
end

for j=1:rob_no %Initial points - robots
    i = i + 1;
    cel = random(i);
    pos_x = mod(cel-1,x);
    pos_y = floor((cel-1)/x);
    initial_point{j}=mean(obs_size*[pos_x pos_x+1 pos_x+1 pos_x; pos_y pos_y pos_y+1 pos_y+1],2);
    if (rob_no > 1)
        final_point{j} = initial_point{j};
    end
end

if (rob_no == 1)
    i = i + 1;
    cel = random(i);
    pos_x = mod(cel-1,x);
    pos_y = floor((cel-1)/x);
    final_point{j}=mean(obs_size*[pos_x pos_x+1 pos_x+1 pos_x; pos_y pos_y pos_y+1 pos_y+1],2);
end

for k=1:obs_no %regions
    i = i + 1;
    cel = random(i);
    pos_x = mod(cel-1,x);
    pos_y = floor((cel-1)/x);
    objects{k}=obs_size*[pos_x pos_x+1 pos_x+1 pos_x; pos_y pos_y pos_y+1 pos_y+1];
end


end




