function Obs = observation_set(OBS_set,N_r,N_p)
%construct possible set of observations for N_r robots
%it is the cartesian product for N_r times of OBS_set of one robot (from which dummy is first removed) 
%in the result, dummy-free proposition will be added as the last one
%N_p - number of propositions (excepting the empty/dummy one - environment)
%N_r - number of robots
%OBS_set - matrix with N_p columns, each row being a possible observation of a robot (some props can overlap, that's why we use this matrix), padded with zeros to end
%last row in OBS_set is for dummy (empty observation), it will be ignored for cartesian product
% Obs will be a matrix with N_p columns

%first, construct all possible observations for N_p props (matrix with 2^N_p lines, last one for dummy) - name this matrix All_Obs
temp_obs=1:N_p;   %atomic propositions, for constructing possible observations (without dummy prop)
temp_cell=mat2cell( repmat(temp_obs,N_p,1) , ones(1,N_p) , length(temp_obs) );  %duplicate observables of transition systems
temp_obs=cartesian_product(temp_cell{:});  %possible observables, on rows (more robots can have the same observables, that's why we use carth product); observables will be labeled with indices of rows (in T.obs)
temp_obs=unique(sort(temp_obs,2),'rows');   %sort rows and remove identical ones (they would have the same observable)
for i=1:size(temp_obs,1)  %modify observables (keep at most one occurence of same prop on each row, and pad with zeros until length 
    obs=unique(temp_obs(i,:));    %keep unique elements on each row, and pad wih zeros
    if length(obs)<N_p %pad with zeros until number of propositions
        obs((end+1):N_p)=0;
    end
    temp_obs(i,:)=obs;
end
temp_obs=unique(temp_obs,'rows');   %again remove identical rows (there may appear new repetitions, due to keeping only unique elements on each row and padding with zeros)
%until now temp_obs has 2^n-1 lines (-1 for the empty set); add a line for the dummy prop (we shouldn't add dummy to other satisfied props, only on a single line)
temp_obs(end+1,:)=[N_p+1 , zeros(1,N_p-1)]; %dummy has index (N_p+1), and pad with zeros after it
All_Obs=temp_obs; %All_Obs contains possible observables of the system (on rows, propositions that can be satisfied; last row for dummy)


%second, find indices in All_Obs of each individual possible observation from OBS_set and construct a row vector (obs_indices) with these indices
obs_indices=zeros(1,size(OBS_set,1));   %number of lines of OBS_set is the number of elements
for i=1:size(OBS_set,1)
    [~,obs_indices(i)]=ismember(OBS_set(i,:),All_Obs,'rows');    %observable
end
obs_indices = setdiff(obs_indices,size(All_Obs,1)); %remove last element, since that corresponds to empty (dummy) observation (N_p+1)
%obs_indices is for one robot (identical robots


%third, take cartesian product of obs_indices for N_r -times
temp_cell=mat2cell( repmat(obs_indices,N_r,1) , ones(1,N_r) , length(obs_indices) );  %duplicate observables for N_r times
temp_obs=cartesian_product(temp_cell{:});  %once having temp_obs, use same syntax as for the first step
temp_obs=unique(sort(temp_obs,2),'rows');   %sort rows and remove identical ones (they would have the same observable)
Obs=zeros(size(temp_obs,1),N_p); %initialize the matrix to be returned (even if it has more rows than necessary)
for i=1:size(temp_obs,1)  %compute and store observables from their indices from each row of temp_obs
    obs_i=unique(temp_obs(i,:));    %keep unique observation indices
    obs=[]; %this will store current observed props for corresponding indices of row i of temp_obs (or obs=zeros(1,0); %to force to be row vector)
    for j=1:length(obs_i)   %add observed props "OBS_set(obs_i(j),:)" to existing set for this row of temp_obs
        obs = union(obs, All_Obs(obs_i(j),:));
    end
    obs=setdiff(obs,0); %obs is already unique and sorted; remove zero elements (we need only propositions)
    obs = reshape(obs,1,[]);
    if length(obs)<N_p %pad with zeros until number of propositions
        obs((end+1):N_p)=0;
    end
    Obs(i,:)=obs;   %add as a row in initialized matrix Obs
end
Obs=unique(Obs,'rows');   %remove identical rows
Obs(end+1,:)=[N_p+1 , zeros(1,N_p-1)]; %dummy has index (N_p+1), and pad with zeros after it
