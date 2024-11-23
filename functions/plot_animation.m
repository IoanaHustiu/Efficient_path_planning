function plot_animation(T,Pre,Post,m,sigma,env_bounds)

RobotInitPlaces = T.R0;
m0 = m(:,1);
mf = m(:,end);
number_intermediary_markings = size(sigma,2); %number of intermediaty markings

     plot_environment(setdiff(find(m0),find(mf)),T.props,T.obstacles,env_bounds,T.Vert); %plot partition
    title(sprintf('Iteration %d',0));
    for i = 1 : number_intermediary_markings
        if (i == 1)
            current_marking = m0;
        else
            current_marking = m(:,i);
            %                             current_marking = round(sol2((nplaces+ntrans)*(i-2)+1: (nplaces+ntrans)*(i-2)+nplaces ));
        end
        plot_environment(setdiff(find(current_marking),find(mf)),T.props,T.obstacles,env_bounds,T.Vert); %partition
        title(sprintf('Iteration %d',i));
        current_sigma = sigma(:,i);
        [feasible_sigma,Rob_places,~,Rob_positions] = sigma2trajectories(Pre,Post,current_marking,current_sigma,RobotInitPlaces);
        if ~feasible_sigma
            error('something wrong happened!');
        end
        RobotInitPlaces = Rob_positions;
        for j = 1 : length(Rob_places)
            traj = Rob_places{j};
            for k = 1 : length(traj)-1
                plot([T.centr{traj(k)}(1) T.centr{traj(k+1)}(1)],[T.centr{traj(k)}(2) T.centr{traj(k+1)}(2)],'-','LineWidth',1,'Color','magenta');
            end
            if (length(traj) > 1)
                init = T.centr{traj(1)};
                fill([init(1)-1.5 init(1)+1.5 init(1)+1.5 init(1)-1.5],[init(2)-1.5 init(2)-1.5 init(2)+1.5 init(2)+1.5],'magenta','EdgeColor','none');
                init = T.centr{traj(end)};
                fill([init(1)-1.5 init(1)+1.5 init(1)+1.5 init(1)-1.5],[init(2)-1.5 init(2)-1.5 init(2)+1.5 init(2)+1.5],'magenta','EdgeColor','none');
            end
            %final destinations already reached
            reached = intersect(find(current_marking),find(mf));
            for k=1:length(reached)%final points reached
                fill(T.Vert{reached(k)}(1,:),T.Vert{reached(k)}(2,:),'green','LineWidth',1);%,'LineStyle','--');,'FaceAlpha',0.5
            end
        end
    end

end