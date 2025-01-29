function T = create_partition_obstacles(regions,o_cells,bounds,varargin)
%partition and observation set; OBS_set is for one robot, i.e. disjoint regions cannot be simmultaneously observed

[C,adj]=grid_decomposition_regions(env_bounds);

T.Q=1:length(C);
T.Vert=C;
T.adj=adj;

for i=1:length(o_cells)
    T.obstacles{i} = T.Vert{o_cells(i)};
end

%find cells corresponding to each defined region (proposition): 
T.props=cell(1,length(regions));  %props{i} will be row vector with indices of cells included in proposition(region) i
for i=1:length(regions)
    ind_obs=find(sum(T.OBS_set==i , 2));  %indices of observables containing prop. i
    T.props{i}=find(ismember(T.obs,ind_obs)); %searched cells
    idx=0;
%     idx = length(find(o_cells<T.props{i})); %find the number of obstacles that are before the current region

    T.props{i} = T.props{i} - idx; %adjust the number of cell according to the environment
end

% T.free_sp = setdiff(1:length(T.Q)-length(o_cells),unique([T.props{:}]));
T.free_sp = setdiff(T.Q,unique([T.props{:}]));    %indices of free cells (not belonging to any region of interest)
% T.free_sp(end-length(o_cells)+1:end)=[];

%initial robot positions are available, find corresponding cells
if nargin==4
    x0=varargin{1};
    N_r = length(x0);
    T.R0=zeros(1,N_r);
    for r=1:N_r
        for i=1:length(T.Q) %find initial discrete state
            if inpolygon(x0{r}(1,1),x0{r}(2,1),T.Vert{i}(1,:),T.Vert{i}(2,:))
%                 idx = length(find(o_cells<i));
%                 T.R0(r)=i-idx;
                T.R0(r) = i;
                break;
            end
        end
    end

end


T.adj(o_cells,:) = 0;
T.adj(:,o_cells) = 0;
T.obs(o_cells) = 0;

for i=1:length(o_cells)
    idx = find(T.free_sp==o_cells(i));
    if ~isempty(idx)
        T.free_sp(idx) = [];
    end
end

% T.Q(end-length(o_cells)+1:end)=[]; 
% T.Vert(o_cells) = [];
% T.adj(o_cells,:) = []; T.adj(:,o_cells) = [];
% T.obs(o_cells) = [];
% T.mid_X(o_cells,:) = []; T.mid_X(:,o_cells) = [];
% T.mid_Y(o_cells,:) = []; T.mid_Y(:,o_cells) = [];