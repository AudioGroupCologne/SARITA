function [nodes, radius] = getEigenmike64Nodes()
% import json file from sparta array2sh ver.1.7.0
% return nodes in radiant, Azimuth and Elevation angles
    fid = fopen('em64.json');
    raw = fread(fid, inf);
    str = char(raw');
    fclose(fid);
    data = jsondecode(str);
    nodes = [[data.GenericLayout.Elements.Azimuth]', ...
             [data.GenericLayout.Elements.Elevation]'];
    
    % make azimuth in range 0 to 360
    nodes(:, 1) = mod(nodes(:, 1) + 360, 360);
    nodes = deg2rad(nodes);
    radius = 0.042;
end