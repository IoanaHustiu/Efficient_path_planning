function env_h = plot_environment(initial_regions,final_regions,obstacles,env_bounds,cells)
%ver. dec.2015

env_h=figure(); %figure handle
axis([0 env_bounds 0 env_bounds]); %axis(world_dim);
hold on
set(gca,'Box','on');
set(gca,'XTick',[],'YTick',[]);
for i=1:length(initial_regions) %initial points
    fill(cells{initial_regions(i)}(1,:),cells{initial_regions(i)}(2,:),'red','LineWidth',1);%,'LineStyle','--');
%    centr=mean(cells{initial_regions(i)},2)';
%    text(centr(1),centr(2),sprintf('y_{%d}',i),'HorizontalAlignment','center','Color','w','FontSize',12,'FontWeight','bold','FontName','Times New Roman');
end

for i=1:length(final_regions)%final points
    fill(cells{final_regions{i}}(1,:),cells{final_regions{i}}(2,:),'blue','LineWidth',1);%,'LineStyle','--');,'FaceAlpha',0.5
%    centr=mean(cells{final_regions{i}},2)';
%    text(centr(1),centr(2),sprintf('o_{%d}',i),'HorizontalAlignment','center','Color','w','FontSize',12,'FontWeight','bold','FontName','Times New Roman');
end

for i=1:length(obstacles) %obstacles
    fill(obstacles{i}(1,:),obstacles{i}(2,:),'black','LineWidth',1);%,'LineStyle','--');'FaceAlpha',0.5
%    centr=mean(cells{obstacles(i)},2)';
%    text(centr(1),centr(2),sprintf('y_{%d}',i),'HorizontalAlignment','center','Color','white','FontSize',12,'FontWeight','bold','FontName','Times New Roman');
end

C=cells;
%represent cells:
for i=1:length(C)
    fill(C{i}(1,:),C{i}(2,:),'w','FaceAlpha',0.5);
end

%write cell number
for i=1:length(C)
    centr=mean(C{i},2)';
%    text(centr(1),centr(2),sprintf('p_{%d}',i),'HorizontalAlignment','center','Color','k','FontSize',5,'FontAngle','italic','FontName','TimesNewRoman');
end
set(gca,'Box','on','XTick',[],'YTick',[]);
end

