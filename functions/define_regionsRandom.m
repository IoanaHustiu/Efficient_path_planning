function [regions,varargout] = define_regionsRandom(env_bounds,reg_number,varargin)
%ver_dec.2015
%env_bounds has the form [x_min,x_max,y_min,y_max]
%reg_number - number of regions to define via mouse clicks
%can read initial position of a given number of robots - if there are three input arguments, last one is the number of robots

init_wrld_h=figure(); %handle al figurii cu lumea initiala
axis(env_bounds);
hold on
grid
regions=cell(1,reg_number);
%read obstacle's vertices
for i=1:reg_number   %i = object number
    for j=1:3
        x = 16*rand;
        y = 16*rand;
        plot(x,y,'.k')
        regions{i}(:,j)=[x;y];
    end

    %creating convex obstacles & drawing them
    k=convhull(regions{i}(1,:),regions{i}(2,:));
    regions{i}=regions{i}(:,k(1:length(k)-1));
    pause(0.3)
    fill(regions{i}(1,:),regions{i}(2,:),'k','FaceAlpha',0.5); %or functia patch (similara cu fill)
end

%if desired, randomly generate initial position of robots
if nargin == 3
    N_r=varargin{1};    %number of robots
    x0=cell(1,N_r);
    for r=1:N_r
        x0{r} = 16*rand(2,1);
        plot(x0{r}(1,1),x0{r}(2,1),'*r');
    end
    varargout{1}=x0;
end

