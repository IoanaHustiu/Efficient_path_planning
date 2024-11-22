function env_h = plot_environment_obstacles(regions,obstacles,env_bounds,reg_edges,varargin)
%ver. dec.2015

env_h=figure(); %figure handle
axis(env_bounds); %axis(world_dim);
hold on
set(gca,'Box','on');
set(gca,'XTick',[],'YTick',[]);
for i=1:length(regions)
    fill(regions{i}(1,:),regions{i}(2,:),'k','FaceAlpha',0.5,'EdgeColor',reg_edges{mod(i-1,length(reg_edges))+1},'LineWidth',3);%,'LineStyle','--');
%     fill(regions{i}(1,:),regions{i}(2,:),'k','FaceAlpha',0.5,'EdgeColor',reg_edges{mod(i-1,length(reg_edges))+1},'LineWidth',2);%,'LineStyle','--');
    % fill(regions{i}(1,:),regions{i}(2,:),'k','FaceAlpha',0.5,'EdgeColor','y','LineWidth',1);%,'LineStyle','--');
    centr=mean(regions{i},2)';
%     if nargin<4
%         text(centr(1),centr(2),sprintf('\\Pi_{%d}',i),'HorizontalAlignment','center','Color','w','FontSize',12,'FontWeight','bold','FontName','Times New Roman');
        text(centr(1),centr(2),sprintf('y_{%d}',i),'HorizontalAlignment','center','Color','w','FontSize',16,'FontWeight','bold','FontName','Times New Roman');
%     end
end

for i=1:length(obstacles)
    fill(obstacles{i}(1,:),obstacles{i}(2,:),'k','FaceAlpha',0.75,'EdgeColor','k','LineWidth',3);%,'LineStyle','--');
%     fill(regions{i}(1,:),regions{i}(2,:),'k','FaceAlpha',0.5,'EdgeColor',reg_edges{mod(i-1,length(reg_edges))+1},'LineWidth',2);%,'LineStyle','--');
    % fill(regions{i}(1,:),regions{i}(2,:),'k','FaceAlpha',0.5,'EdgeColor','y','LineWidth',1);%,'LineStyle','--');
    centr=mean(obstacles{i},2)';
%     if nargin<4
%         text(centr(1),centr(2),sprintf('\\Pi_{%d}',i),'HorizontalAlignment','center','Color','w','FontSize',12,'FontWeight','bold','FontName','Times New Roman');
        text(centr(1),centr(2),sprintf('o_{%d}',i),'HorizontalAlignment','center','Color','w','FontSize',16,'FontWeight','bold','FontName','Times New Roman');
%     end
end

switch nargin
    case 4
        %do nothing (obstacles already plotted)
    case 5  %argument C (cells) - plot cell decomposition
        C=varargin{1};
        %represent cells:
        for i=1:length(C)
            fill(C{i}(1,:),C{i}(2,:),'w','FaceAlpha',0.5);
        end
        
        %write cell number
        for i=1:length(C)
            centr=mean(C{i},2)';
            text(centr(1),centr(2),sprintf('p_{%d}',i),'HorizontalAlignment','center','Color','k','FontSize',5,'FontAngle','italic','FontName','TimesNewRoman');
        end
        set(gca,'Box','on','XTick',[],'YTick',[]);
end
end

