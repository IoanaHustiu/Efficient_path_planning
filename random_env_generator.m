% Function to create Random Enviroment
function [objects,obstacles,o_cells,initial_point,final_point,adj,rem_cells] = random_env_generator(obs_no,N_o,obs_size,rob_no,limits)

objects={};
objects_bound={};

obstacles={};
obstacles_bound={};

initial_point={};
initial_bound={};

final_point={};
final_bound={};

%% initial positions for the robots + ROI

x = (limits(2)-limits(1)) / obs_size;
y = (limits(4)-limits(3)) / obs_size;

random = [];
o_cells = [];
obstac = [];

% generate obstacles in the environment
flag = 1;

while flag
    while length(random) < N_o
        obstac = unique(randi([1 x*y], 1, N_o-length(random)),'stable');
        random = unique([random obstac],'stable');
    end
    
% create graph of the remaining cells
    C = cell(1,x*y);
    adj = sparse(eye(x*y)); %self transitions in every cell
    
    for j=0:y-1
        for i=0:x-1
            C{j*x+i+1} = obs_size*[i i+1 i+1 i; j j j+1 j+1];
        end
    end

nr_reg = length(C);

    for i=1:nr_reg
        for j=i+1:nr_reg %go only through higher indices (adjacency matrix is symmetric)
%             if ~any(obstac==i) && ~any(obstac==j)
                if (norm(mean(C{i},2)-mean(C{j},2)) <= 10^10*eps+obs_size) %two common vertices means adjacency
                    adj(i,j) = 1;
                    adj(j,i) = 1; %symmetric
                end
%             end
        end
    end

    adj(random,:) = 0;
    adj(:,random) = 0;

% find connected components
    G = graph(adj);
    [str_bins, binsize] = conncomp(G);
    I = binsize(str_bins) == max(binsize);

    rem_cells = find(I==1);

    if floor(max(binsize)/2) < obs_no
        fprintf('Not enough connected cells. New environment is generated.\n')
    else
        flag = 0;
    end
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

    %draw the region
    % % %     pause(0.001);
    % % %     fill(objects{i}(1,:),objects{i}(2,:),'b-','FaceAlpha',0.5); %or functia patch (similara cu fill)
end

%% Plotting points
% % %
% % % for i=1:length(initial_point)
% % %     plot(initial_point{i}(1),initial_point{i}(2),'color',colors{i},'marker','o','Linewidth',2);
% % %     text((initial_point{i}(1)+0.2),(initial_point{i}(2)+0.2),{num2str(i)});
% % %     hold on;
% % % end
% % %
% % % if (rob_no == 1)
% % %     for i=1:length(final_point)
% % %         plot(final_point{i}(1),final_point{i}(2),'color','k','marker','x','Linewidth',2);
% % %         text((final_point{i}(1)+0.2),(final_point{i}(2)+0.2),{char(64+i)});
% % %         hold on;
% % %     end
% % % end
% % % pause(0.1)

end




