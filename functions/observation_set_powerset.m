function Obs = observation_set_powerset(N_p)
%construct possible set of observations (dummy-free proposition will be the last one)
%N_p - number of propositions (excepting the dummy one-environment) that can be simultaneously satisfied 
%consider maximum number of possible simultaneous observations (assume more robots or/and that all props. can overlap)
% Obs will correspond to the power set of 1:N_p
% Obs will be a matrix with 2^N_p rows and N_p columns

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
Obs=temp_obs; %Obs contains possible observables of the system (on rows, propositions that can be satisfied; last row for dummy)
