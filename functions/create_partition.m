function T = create_partition(regions,bounds,varargin)
%ver. dec.2015
%partition and observation set; OBS_set is for one robot, i.e. disjoint regions cannot be simmultaneously observed

% [C,adj,OBS_set,obs,mid_X,mid_Y]=triangular_decomposition_regions(regions,bounds);
[C,adj,OBS_set,obs,mid_X,mid_Y] = rmt_grid_decomposition_regions(regions,bounds);

T.Q=1:length(C);
T.Vert=C;
T.adj=adj;
T.OBS_set=OBS_set;
T.obs=obs;
T.mid_X=mid_X;
T.mid_Y=mid_Y;

%find cells corresponding to each defined region (proposition): 
T.props=cell(1,length(regions));  %props{i} will be row vector with indices of cells included in proposition(region) i
for i=1:length(regions)
    ind_obs=find(sum(T.OBS_set==i , 2));  %indices of observables containing prop. i
    T.props{i}=find(ismember(T.obs,ind_obs)); %searched cells
end

T.free_sp = setdiff(T.Q , unique([T.props{:}]));    %indices of free cells (not beonging to any region of interest)

%initial robot positions are available, find corresponding cells
if nargin==3
    x0=varargin{1};
    N_r = length(x0);
    T.R0=zeros(1,N_r);
    for r=1:N_r
        for i=1:length(T.Q) %find initial discrete state
            if inpolygon(x0{r}(1,1),x0{r}(2,1),T.Vert{i}(1,:),T.Vert{i}(2,:))
                T.R0(r)=i;
                break;
            end
        end
    end

end
