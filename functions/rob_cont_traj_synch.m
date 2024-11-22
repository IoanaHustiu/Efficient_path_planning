function [rob_traj, synch_points] = rob_cont_traj_synch(T,Runs,Synch,x0)

rob_traj=cell(1,length(Runs));
for r=1:length(Runs)    %for each robot
    rob_traj{r}=zeros(2,length(Runs{r})+1);
    rob_traj{r}(:,1) = x0{r};
    
    for i=2:length(Runs{r})
        if Runs{r}(i-1)~=Runs{r}(i)   %this should be always true
            rob_traj{r}(:,i) = [T.mid_X(Runs{r}(i-1),Runs{r}(i)) ; T.mid_Y(Runs{r}(i-1),Runs{r}(i))];  %middle point of segment shared by two successive states
        end
    end
    
    if length(Runs{r})>1    %robot r has moved
        rob_traj{r}(:,end) = mean(T.Vert{Runs{r}(end)},2);   %centroid of last state
    else
        rob_traj{r}(:,end) = x0{r};    %robot stays in initial point
    end
end

synch_points=cell(1,length(Synch));
for r=1:length(Runs)    %for each robot
    synch_points{r}=zeros(2,length(Synch{r}));
    for i=1:length(Synch{r})
        c_ind=Synch{r}(i); %synchronize with others when entering the "Synch{r}(i)"-th cell from run (if robot moves; otherwise, synchronize in current position)
        if c_ind==1 %first cell -> initial position
            synch_points{r}(:,i) = x0{r};
        elseif Runs{r}(c_ind-1)~=Runs{r}(c_ind) %robot moves; synchronize when entering
            synch_points{r}(:,i) = [T.mid_X(Runs{r}(c_ind-1),Runs{r}(c_ind)) ; T.mid_Y(Runs{r}(c_ind-1),Runs{r}(c_ind))];  %middle point of segment
        elseif Runs{r}(c_ind-1)==Runs{r}(c_ind) %robot doesn't move now; synchronize in previous point
            synch_points{r}(:,i) = synch_points{r}(:,i-1); %or can choose other point inside current cells
        end
    end
end
