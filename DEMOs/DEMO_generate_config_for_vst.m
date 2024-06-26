% dependencies:
% aktools, sfs, sofia, sfe 
% generate configs for VST plugin
clear; 
close all; 
clc;

addpath(genpath('../dependencies'));
addpath('../SARITA_code');

path_to_configs = 'sarita_vst_configs';

N_ups = 20;
target_grid = sph_grids.get_sampling_grid('lebedev', N_ups);
src_grids = {'zylia', 'eigenmike32', 'eigenmike64', 'hosma'};

for g = 1:numel(src_grids)
    [source_grid, ~, N, radius] = sph_grids.get_sampling_grid(src_grids(g));

    config.N = N;
    config.N_ups = N_ups;
    config.radius = radius;
    config.source_grid_name = string(src_grids(g));
    config.path = path_to_configs; 
    if ~exist(config.path, 'dir') 
        mkdir(config.path); 
    end
    array_input = zeros(size(source_grid, 1), 2048);
    [~] = Sarita_upsampling_faster(array_input, source_grid(:, 1:2), target_grid(:, 1:2), radius, ...
                                 'frame_length', 4096, 'frame_overlap', 128, ...
                                 'gen_vst_config', 1, 'config', config);
end

