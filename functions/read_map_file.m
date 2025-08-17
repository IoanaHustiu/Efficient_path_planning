function grid = read_map_file(filename)
    fid = fopen(filename, 'r');
    line = fgetl(fid);
    
    while ischar(line)
        if startsWith(line, 'height')
            height = str2double(extractAfter(line, 'height '));
        elseif startsWith(line, 'width')
            width = str2double(extractAfter(line, 'width '));
        elseif strcmp(line, 'map')
            break;
        end
        line = fgetl(fid);
    end

    grid = zeros(height, width);  % initialize
    for i = 1:height
        line = fgetl(fid);
        for j = 1:width
            if line(j) == '.'
                grid(i,j) = 0; % free
            else
                grid(i,j) = 1; % obstacle
            end
        end
    end

    fclose(fid);
end