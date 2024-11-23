function rob_traj = compute_trajectories(T,Pre,Post,m,sigma,RobotInitPositions)

N_r = length(T.R0);
rob_traj = RobotInitPositions; %initialize the trajectories with initial positions of the robots
number_intermediary_markings = size(sigma,2); %number of intermediaty markings
RobotInitPlaces = T.R0; %initialize the initial places of the robots

for i=1:number_intermediary_markings
    [feasible_sigma,Rob_places,~,Rob_positions] = sigma2trajectories(Pre,Post,m(:,i),sigma(:,i),RobotInitPlaces); % for all the sigma vectors that are feasible, choose one for which the number of transitions is minimum
    if feasible_sigma
        for r=1:N_r
            for i1=1:size(Rob_places{r},2)-1
                rob_traj{r}(:,end+1) = [T.mid_X(Rob_places{r}(i1),Rob_places{r}(i1+1)); T.mid_Y(Rob_places{r}(i1),Rob_places{r}(i1+1))];
            end
            if size(Rob_places{r},2)>1
                rob_traj{r}(:,end+1) = mean(T.Vert{Rob_places{r}(end)},2); %center of the final position
            end
        end
        RobotInitPlaces= Rob_positions;
    end
end

end