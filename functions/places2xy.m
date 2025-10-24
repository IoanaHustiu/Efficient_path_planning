function XY = places2xy(places, T)
% PLACES2XY  Converts a vector of place indices [p1 p2 ...]
%             into an N×2 matrix of coordinates [x y].
%
%   XY = places2xy(places, T)
%
%   Each T.centr{i} must be a 1×2 vector [x y].

    places = places(:);
    n = numel(places);
    XY = zeros(n, 2);

    if isfield(T,'centr') && numel(T.centr) >= max(places)
        for k = 1:n
            XY(k, :) = T.centr{places(k)};
        end
    else
        [H, W] = size(T.map2D);
        lin = T.rem_cells(places);
        [r, c] = ind2sub([H, W], lin);
        XY = [c, r];
    end
end

