function [td_data, grid_data, fs, radius, td_struct, N_grid] = sfe_import_sofa(filename, find_grid_type)
if nargin < 2
    find_grid_type = 0;
end
try
    sfob = SOFAload(filename);
catch
    SOFAstart();
    sfob = SOFAload(filename);
end

if strcmp(sfob.GLOBAL_SOFAConventions, 'SingleRoomDRIR')
    fprintf('sfe_import_sofa: Load SingleRoomDRIR %s\n', string(filename))
    td_data = squeeze(sfob.Data.IR);
    grid_data = sfob.ReceiverPosition;
    radius = grid_data(1, 3);
    if sfob.ReceiverPosition_Type == 'cartesian'
        [az, el, r] = cart2sph(grid_data(:, 1), grid_data(:, 2), grid_data(:, 3));
        az = mod(az, pi*2);
        col = pi/2-el;
        grid_data = [az, col];
        radius = r(1);
    elseif sfob.ReceiverPosition_Type == 'spherical'
        grid_data = deg2rad(grid_data);
        grid_data(:, 2) = pi/2 - grid_data(:, 2);
    end
    grid_data = grid_data(:, 1:2);
    
elseif strcmp(sfob.GLOBAL_SOFAConventions, 'MultiSpeakerBRIR')
    fprintf('sfe_import_sofa: Load MultiSpeakerBRIR %s\n', string(filename))
    td_data = permute(squeeze(sfob.Data.IR), [1, 3, 2]);
    grid_data = deg2rad(sfob.ListenerView(:, 1:2));
    radius = sfob.ListenerView(1, 3);
    
elseif strcmp(sfob.GLOBAL_SOFAConventions, 'SimpleFreeFieldHRIR')
    fprintf('sfe_import_sofa: Load SimpleFreeFieldHRIR %s\n', string(filename))
    td_data = permute(squeeze(sfob.Data.IR), [1, 3, 2]);
    grid_data = sfob.SourcePosition;
    grid_data(:, 2) = 90-grid_data(:, 2);
    grid_data(:, 1:2) = deg2rad(grid_data(:, 1:2));
    radius = grid_data(1, 3);
    grid_data = grid_data(:, 1:2);
end

if find_grid_type
    [N_grid, weights, name] = sph_grids.find_grid_type(grid_data);
    if ~isempty(N_grid)
        warning(sprintf('Found grid: %s N=%d, --> adopt parameters', name, N_grid));
        grid_data(:, 3) = weights;
    end
else
    N_grid = [];
end

fs = sfob.Data.SamplingRate;
disp('sfe_import_sofa: No grid weights stored in SOFA files so far')

% set up sofia time domain struct
td_struct.impulseResponses = td_data;
td_struct.FS               = fs;
td_struct.averageAirTemp   = 20;
td_struct.radius           = radius;
end

