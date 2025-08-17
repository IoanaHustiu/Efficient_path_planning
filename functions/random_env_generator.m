function obstacles = random_env_generator(N_o,N_regions)
% Function to create Random Enviroment


obstacles = [];
% generate obstacles in the environment

while (length(obstacles) < N_o)
    obstac = unique(randi([1 N_regions], 1, N_o-length(obstacles)),'stable');
    obstacles = unique([obstacles obstac],'stable');
end

return

end




