function plot_trajectories(rob_traj,rob_color)

N_r = size(rob_traj,2);

% for r=1:N_r    %plot initial positions of robots
%     plot(x0{r}(1,1),x0{r}(2,1),rob_color{mod(r-1,length(rob_color))+1},'Marker','o','LineWidth',1.5);
% end

for r=1:N_r    %plot trajectories of robots
    plot(rob_traj{r}(1,1),rob_traj{r}(2,1),rob_color{mod(r-1,length(rob_color))+1},'Marker','o','LineWidth',1.5);
    plot(rob_traj{r}(1,1:end),rob_traj{r}(2,1:end),rob_color{mod(r-1,length(rob_color))+1},'LineWidth',1.5);
    plot(rob_traj{r}(1,end),rob_traj{r}(2,end),rob_color{mod(r-1,length(rob_color))+1},'Marker','x','LineWidth',1.5);
end

end