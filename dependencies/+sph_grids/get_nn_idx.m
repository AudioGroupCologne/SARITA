function idx = get_nn_idx(grid_data, az, col)
%
% get the index of sampling point closest to desired azimuth and colatitude
%
% Parameter:
% ------------
%
%   grid_data: azimuth and coltitude angles of the sperical grid
%   az:        azimuth angle in RAD 
%   col:       colatitude angle in RAD
%
% Returns:
% ------------
% 
% idx: index in grid data, closest to az and col
%
%
% (C) 01/2020 Tim Luebeck
    grid_data(:, 1) = mod(grid_data(:, 1)+2*pi, 2*pi);
    grid_data(:, 2) = (pi/2)-grid_data(:, 2);
    
    az = mod(az+2*pi, 2*pi);
    col = (pi/2)-col;
    
    [grid_cart(:, 1), grid_cart(:, 2), grid_cart(:, 3)] = ...
            sph2cart(grid_data(:, 1), grid_data(:, 2), ones(size(grid_data, 1), 1));
    [x, y, z] = sph2cart(az, col, 1);
    
    for i = 1:size(x, 1)
        dists(:, i) = sqrt(sum([(grid_cart(:, 1)-x(i)).^2, (grid_cart(:, 2)-y(i)).^2, (grid_cart(:, 3)-z(i)).^2], 2));
    end
    [~, idx] = min(dists);
end
