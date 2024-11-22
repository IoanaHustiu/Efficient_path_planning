function [regions,varargout] = define_regions(env_bounds,reg_number,varargin)
%ver_dec.2015
%env_bounds has the form [x_min,x_max,y_min,y_max]
%reg_number - number of regions to define via mouse clicks
%can read initial position of a given number of robots - if there are three input arguments, last one is the number of robots

% disp('Definirea unor regiuni convexe (minim 3 varfuri) folosind mouse-ul');

fprintf('\nFor defining a region:\n\t - left-click = pick a vertex\n\t - right-click = pick last vertex\nRegions should be convex and non-overlapping\n\nPress any key to begin.\n')
pause();

% scrsz = get(0,'ScreenSize'); %rezolutia display-ului (pt a maximiza figurile)
% scrsz(4)=scrsz(4)-72;

init_wrld_h=figure(); %handle al figurii cu lumea initiala
% set(init_wrld_h,'Position',scrsz); %maximizam figura
axis(env_bounds);
% title('Obstacle vertices')
hold on
grid
regions=cell(1,reg_number);
%read obstacle's vertices
for i=1:reg_number   %i = object number
    j=1; %j = no. of vertexes for current object
    but=1;
    while but==1
        [x,y,but]=ginput(1);
%                 x=round(x*2)/2;
%                 y=round(y*2)/2;
        plot(x,y,'.k')
        regions{i}(:,j)=[x;y];
        j=j+1;
    end
    
    %creating convex obstacles & drawing them
    k=convhull(regions{i}(1,:),regions{i}(2,:));
    regions{i}=regions{i}(:,k(1:length(k)-1));
    pause(0.3)
    fill(regions{i}(1,:),regions{i}(2,:),'k','FaceAlpha',0.5); %or functia patch (similara cu fill)
end

%if desired, read initial position of robots
if nargin == 3
    N_r=varargin{1};    %number of robots
    x0=cell(1,N_r);
    for r=1:N_r
        x0{r}=zeros(2,1);
        title(sprintf('Click on initial position of robot %d',r));
        [x,y]=ginput(1);
        x0{r}(1,1)=x;   %initial continuous position
        x0{r}(2,1)=y;
        plot(x0{r}(1,1),x0{r}(2,1),'*r');
    end
    varargout{1}=x0;
end

% close(init_wrld_h); %close figure (an environment one will be drawn outside this function)
