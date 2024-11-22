function env_h = plot_environment(regions,env_bounds,reg_edges,varargin)
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

switch nargin
    case 3
        %do nothing (obstacles already plotted)
    case 4  %argument C (cells) - plot cell decomposition
        C=varargin{1};
        %represent cells:
        for i=1:length(C)
            fill(C{i}(1,:),C{i}(2,:),'w','FaceAlpha',0.5);
        end
        
        %write cell number
        for i=1:length(C)
            centr=mean(C{i},2)';
            text(centr(1),centr(2),sprintf('p_{%d}',i),'HorizontalAlignment','center','Color','k','FontSize',10,'FontAngle','italic','FontName','TimesNewRoman');
        end
        set(gca,'Box','on','XTick',[],'YTick',[]);
        
    case 5 %arguments C (cells) and adjacency - plot adjacency graph and cell decomposition
        C=varargin{1};
        adj=varargin{2};
        %represent cells:
        for i=1:length(C)
            fill(C{i}(1,:),C{i}(2,:),'w','FaceAlpha',0.5);
        end

        centr=zeros(length(C),2);   %store centroids
        for i=1:length(C)
            centr(i,:)=mean(C{i},2)';
        end
        gplot(adj,centr,':b');  %represent adjacency graph

        %write cell number
        for i=1:length(C)
            text(centr(i,1),centr(i,2),sprintf('p_{%d}',i),'HorizontalAlignment','center','Color','k','FontSize',10,'FontAngle','italic','FontName','TimesNewRoman');
        end
        set(gca,'Box','on','XTick',[],'YTick',[]);
        
    case 7  %arguments C,adj,middle_X,middle_Y - do not represent adjacency graph
        C=varargin{1};
        adj=varargin{2};
        middle_X=varargin{3};
        middle_Y=varargin{4};
        %represent cells:
        for i=1:length(C)
            fill(C{i}(1,:),C{i}(2,:),'w','FaceAlpha',0.5);
        end

        centr=zeros(length(C),2);   %store centroids
        for i=1:length(C)
            centr(i,:)=mean(C{i},2)';
        end
%         gplot(adj,centr,':b');  %represent adjacency graph

        for i=1:length(C)
            for j=setdiff(find(adj(i,:)),1:i)
                plot(middle_X(i,j),middle_Y(i,j),'*b')
            end
        end
        %write cell number
        for i=1:length(C)
            text(centr(i,1),centr(i,2),sprintf('p_{%d}',i),'HorizontalAlignment','center','Color','k','FontSize',10,'FontAngle','italic','FontName','TimesNewRoman');
        end
        set(gca,'Box','on','XTick',[],'YTick',[]);
        
    otherwise
        fprintf('\nCheck # of arguments for the plot\_environment function\n')
end

