%% SARITA - Spherical Array Interpolation by Time Alignment
%
% function upsampledArrayDataset = Sarita_upsample(sparseArrayDataset, denseSamplingGrid,fTransit,frameLength)
%
% This function performs the spatial upsampling of a sparse spherical array dataset
% to a predefined denseSamplingGrid. 
% 
% The weights for the different directions are determined using Voronioi diagrams 
% A specific method of temporally aligning the impulse responses is
% performed. Further information in XXXXX
%
% This version is much simpler than the normal SARITA algorithm.
% - No subsample accuracy for the Peak detection and shifts
% - No separated low-frequency processing
%
% Parameters:
% ------------
%   sparseArrayDataset:  struct according to sarita_gen_array_struct with 
%                        the following fields:
%                        .P: frequency domain data of the array [M X N],
%                        whith M: number of measurements, N: number of
%                        samples
%                        .samplingGrid: in [M X 2] in DEG and COL, where M
%                                       with M sampling points
%                        .f frequency vector holding the frequency
%                           coresponding to the bins in .P
%                        .FFToversize: oversize
%                        .radius: radius of the array the sound pressure
%                                 has been sampled on
%
%   denseSamplingGrid:   in [M X 2] in DEG and COL, where M
%                                       with M sampling points
%

%
%   frameLength: {default: 32}
%   
%   frameOverlap: (default: frameLength/2);
% 
%   c: {default: 340}
%
% Returns:
% -----------
%
%   upsampledArrayDataset: onesided frequency domain bins
%
% Dependencies:
% -------------
%
% SFS Toolbox >=2.5.0 https://sfs-matlab.readthedocs.io/en/2.5.0/#              
%
% References:
% -------------
%
% [1] 
%
%  (C)   2020 CP,  Christoph P?rschmann
%        TH K?ln - University of Applied Sciences
%        Institute of Communications Engineering
%        Department of Acoustics and Audio Signal Processing
%
% latest update 20.07.23
function drirs_upsampled = Sarita_upsampling_faster(drirs, source_grid, target_grid, radius, varargin)
% check inputs
if nargin < 4
    error('Not enough input arguments');
end
input_params = {'frame_length', 'frame_overlap', 'c', 'fs', 'gen_vst_config', 'config'}; 

default.frame_length   = 32; 
default.frame_overlap  = 0;
default.c              = 340;
default.fs             = 48e3;
default.gen_vst_config = false; % if set, default.config is mandatory!
default.config         = [];

if mod(size(varargin, 2), 2)
    error('Check passed parameter value pairs.')
end

% check inputs
for p_IDX = 1 : 2 : size(varargin, 2)-1
    is_valid = 0;
    for d_IDX = 1:size(input_params, 2)
        if strcmp(input_params{d_IDX}, varargin{p_IDX})
            [~] = evalc([input_params{d_IDX}, '=', 'varargin{p_IDX+1}']);
            is_valid = 1;
        end
    end
    if ~is_valid
        warning('%s is not valid parameter. It will be ignored ', varargin{p_IDX});
    end
end 

% fill all params which were not passed with default values 
for d_IDX = 1:size(input_params, 2)
    value = ['default.' input_params{d_IDX}];
    if ~exist(input_params{d_IDX}, 'var')
        [~] = evalc([input_params{d_IDX}, '=', value]); 
    end
end

frame_length = floor(frame_length/2)*2; % force frameLength to be even
if frame_overlap == 0 || frame_overlap > frame_length/2
    frame_overlap = frame_length/2;
end
clear default varargin

%% end checking arguments 
disp('SARITA -- perform (fast) upsampling of array data')
hannsize = frame_overlap*2;
hannwin = hann(hannsize).';
win = [hannwin(1 : length(hannwin)/2), ...
       ones(1, frame_length-length(hannwin)), ...
       hannwin(length(hannwin)/2+1 : length(hannwin))];

[sparse_grid_cart(:, 1), sparse_grid_cart(:, 2), sparse_grid_cart(:, 3)] = ... 
    sph2cart(source_grid(:, 1), ...
             pi/2 - source_grid(:, 2), ones(length(source_grid(:, 2)), 1));    

 
% allocate some memory

num_neighbors_dense_grid       = zeros(length(target_grid));    % Number of nearest neighbors for each sampling point
idx_neighbors_dense_grid       = zeros(length(target_grid), 4); % Indices of neighbors of each sampling point (max 4)
maxShift_dense                 = zeros(length(target_grid), 1); % Maximal Shift in each sampling point 
weights_neighbors_dense_grid   = zeros(length(target_grid), 4); % Weights of neighbors of each sampling point 
neighbors_combinations         = zeros(2, 0);                   % Array containing all combinations of nearest neighbors
combination_ptr                = zeros(2, 0);                   % Array describing which neighborsCombination is required for each the cross correlations   
 
for dense_idx = 1:size(target_grid, 1)         
    % find nearest neighbors and calculate interpolation weights
    [x, y, z] = sph2cart(target_grid(dense_idx, 1), ...
                         pi/2 - target_grid(dense_idx, 2), ...
                         1);
                     
    [neighborsIndex, weights] = findvoronoi(sparse_grid_cart, [x, y, z]);
  
    % store them
    num_neighbors_dense_grid(dense_idx) = length(neighborsIndex);
    idx_neighbors_dense_grid(dense_idx, 1:num_neighbors_dense_grid(dense_idx)) = neighborsIndex;
    weights_neighbors_dense_grid(dense_idx, 1:num_neighbors_dense_grid(dense_idx)) = weights;
       
    % calculate the maximum timeshifts between all next neighbor
    % candidates and all combinations next neighbors
    for nodeIndex = 2:length(neighborsIndex)   
        angle = acos(dot(sparse_grid_cart(neighborsIndex(1),:), sparse_grid_cart(neighborsIndex(nodeIndex),:))); % get angle between sampling points
        maxShift_dense(dense_idx,nodeIndex-1) = abs(ceil(angle * radius*fs/c)); % calculate the maximum shift                        

        % Determine combinations of nearest neighbors
        compareEntries = ismember(neighbors_combinations, [neighborsIndex(1) neighborsIndex(nodeIndex)]');
        if ~isempty(compareEntries)            
            if (max(compareEntries(1, :) .* compareEntries(2, :))  == 0)                 
                % entry does not exist
                neighbors_combinations = [neighbors_combinations, [neighborsIndex(1) neighborsIndex(nodeIndex)]'];
                combination_ptr = [combination_ptr, [length(neighbors_combinations(1,:)) 1]'];
            else
                % entry already exist, just store the correct reference in
                % combination_ptr
               [~, position]=max(compareEntries(1, :).*compareEntries(2, :));
                if neighbors_combinations(:, position) == [neighborsIndex(1), neighborsIndex(nodeIndex)]' 
                    % entry exists in the same order 
                    combination_ptr = [combination_ptr, [position 1]'];
                else
                    % entry exists in the inverted order
                    combination_ptr = [combination_ptr, [position -1]'];
                end
            end
        else
            % no element in list yet. First element is put in list
            neighbors_combinations = [neighbors_combinations [neighborsIndex(1) neighborsIndex(nodeIndex)]'];
            combination_ptr = [combination_ptr [length(neighbors_combinations(1,:)) 1]'];
        end
        
    end
end      
maxShiftOverall = max(max(abs(maxShift_dense))); % get the maximal possible time shift, required to increase buffer at the end

if gen_vst_config
    % check for necessary parameters
    assert(isfield(config, 'N'));
    assert(isfield(config, 'source_grid_name'));
    
    % setup struct of config parameters
    config.fs = fs;
    config.frame_length = frame_length;
    config.frame_overlap = frame_overlap;
    config.target_grid = target_grid;
    config.maxShiftOverall = maxShiftOverall;
    config.neighbors_combinations = neighbors_combinations;
    config.idx_neighbors_dense_grid = idx_neighbors_dense_grid;
    config.num_neighbors_dense_grid = num_neighbors_dense_grid;
    config.idx_neighbors_dense_grid = idx_neighbors_dense_grid;
    config.weights_neighbors_dense_grid = weights_neighbors_dense_grid;
    config.maxShift_dense  = maxShift_dense;
    config.combination_ptr = combination_ptr;
    
    % drop config to binary file
    gen_sarita_vst_config(config);
end

drirs = [drirs, zeros(size(drirs, 1), maxShiftOverall)]; 

% allocate memory for upsampled drirs
drirs_upsampled = zeros(length(target_grid), length(drirs(1, :)) + maxShiftOverall*2);

number_of_frames = floor(length(drirs(neighborsIndex, :))/(frame_length-frame_overlap))-1;
for frameCounter = 1:number_of_frames % Loop over all frames                       
    startTab = (frameCounter-1)*(frame_length-frame_overlap)+1;
    endTab = startTab+frame_length-1;
    irsFrame = repmat(win, [size(drirs, 1), 1]) .* drirs(:, startTab:endTab); % make compatible with older Matlabversions with repmat
    
    % In each frame the cross-correlation required for the upasmpling are
    % determined 
    % In case of a 350er Grid about 20 % of comp.power
    for crossCorrelationIndex = 1:length(neighbors_combinations)  
         frameOne = irsFrame(neighbors_combinations(1, crossCorrelationIndex), :);
         frameTwo = irsFrame(neighbors_combinations(2, crossCorrelationIndex), :);
         correlationsFrame(:, crossCorrelationIndex) = xcorr(frameOne, frameTwo); 
    end
    neighborsIndexCounter = 0; %Counter which entry in combination_ptr is to be assessed 
    for dirIndex = 1:size(target_grid, 1)    
        timeShiftMean = 0;
        currentTimeShift = 0;
        % Get nearst neighbors, weights and maxShift for actual direction, can
        % later be directly addressed in following lines if desired
        neighborsIndex = idx_neighbors_dense_grid(dirIndex,1:num_neighbors_dense_grid(dirIndex));
        weights = weights_neighbors_dense_grid(dirIndex,1:num_neighbors_dense_grid(dirIndex));
        maxShift = maxShift_dense(dirIndex,1:num_neighbors_dense_grid(dirIndex)-1);            
        neighborsIRs = irsFrame(neighborsIndex,:); % get all next neighbour irs of one frame and perform windowing              
          
        for nodeIndex = 2:length(neighborsIndex)            
            neighborsIndexCounter=neighborsIndexCounter+1;
            correlation = correlationsFrame(:,combination_ptr(1,neighborsIndexCounter));           
            if combination_ptr(2,neighborsIndexCounter)==-1
            	correlation = correlation(end:-1:1);
            end   
            
            % look for maximal value in the crosscorrelated IRs only in the relevant area                            
            correlation = correlation(frame_length - maxShift(nodeIndex-1):frame_length + maxShift(nodeIndex-1));
            [~, maxpos] = max(correlation);
            currentTimeShift(nodeIndex) = (maxpos-(length(correlation)+1)/2);            
            timeShiftMean = timeShiftMean + currentTimeShift(nodeIndex) * weights(nodeIndex);
        end
        
        % align every block according to the calculated time shift, weight
        % and sum up
        % The following loop needs about 40 % comp.power
        for nodeIndex = 1:length(neighborsIndex)
            currentBlock = neighborsIRs(nodeIndex, :) * weights(nodeIndex);
            timeShiftFinal = round(-timeShiftMean + currentTimeShift(nodeIndex) + maxShiftOverall);    % As maxShiftOverall is added, timeShiftFinal will always be positive        
            if timeShiftFinal < 0 % Added 22.12.2021 to assure that timeShiftFinal does not become negative
                timeShiftFinal = 0;
            end
            drirs_upsampled(dirIndex, startTab + timeShiftFinal:endTab + timeShiftFinal) = ...
                drirs_upsampled(dirIndex, startTab + timeShiftFinal:endTab + timeShiftFinal) + currentBlock;
        end
    end
end

drirs_upsampled = drirs_upsampled(:, maxShiftOverall + 1 : length(drirs_upsampled(1, :)) - 2 * maxShiftOverall);
disp('... done')
end 

function gen_sarita_vst_config(config)
    % Method to drop sarita parameters for rendering in VST plugin
    %
    % Gary Grutzek & Tim L?beck / Feb. 2023
    % 
    % /*
    %  * config data struct
    %  */
    % struct SaritaConfig {
    %     // header
    %     uint32_t fs;
    %     uint32_t N;                 // order of source grid
    %     uint32_t NUpsampling;       // order of target grid
    %     uint32_t NRendering;        // rendering order
    %     float radius;               // radius of array
    %     uint32_t denseGridSize;
    %     uint32_t maxShiftOverall;
    %     uint32_t neighborCombLength;   // neighborCombinations size = neighborCombLength * 2
    %     uint32_t idxNeighborsDenseLen; // idxNeighborsDense array size = idxNeighborsDenseLen * dense grid size
    %     uint32_t combinationsPtrLen;   // length of combinations pointer, y is always 2
    % 
    %     // data
    %     uint8_t** neighborCombinations; // Array containing all combinations of nearest neighbors
    %     uint8_t* numNeighborsDense;     // Number of nearest neighbors for each sampling point
    %     uint8_t** idxNeighborsDense;    // Indices of neighbors of each sampling point
    %     float** weightsNeighborsDense;  // Weights of neighbors of each sampling point
    %     uint8_t** maxShiftDense;
    %     int8_t** combinationsPtr;       // Array describing which neighbors combination is required for each cross correlations
    %     float** denseGrid;              // Az, El and weight of each target sensor
    % } cfg;

    if ~isfield(config, 'path') || isempty(config.path)
        fid = fopen(sprintf('Sarita_%s_N%d.cfg', ...
                        config.source_grid_name, config.N), 'wb');
    else
        fid = fopen(sprintf('%s/Sarita_%s_N%d.cfg', ...
                            config.path, ...
                            config.source_grid_name, config.N), 'wb');
    end
    
    % write header
    fwrite(fid, config.fs, 'uint');
    fwrite(fid, config.N, 'uint');                                 % order of source grid
    fwrite(fid, config.N_ups, 'uint');                             % order of target grid
    fwrite(fid, config.radius, 'float')                            % radius of array
    fwrite(fid, size(config.target_grid, 1), 'uint');              % target grid size
    fwrite(fid, config.maxShiftOverall, 'uint');                   % max shift in samples
    fwrite(fid, max(size(config.neighbors_combinations)), 'uint'); % neighborCombinations size
    fwrite(fid, size(config.idx_neighbors_dense_grid, 2), 'uint'); % idxNeighborsDense array size
    fwrite(fid, size(config.combination_ptr, 2), 'uint');          % length of combinations pointer

    % payload
    fwrite(fid, config.neighbors_combinations, 'uint8');
    fwrite(fid, config.num_neighbors_dense_grid(:,1), 'uint8');
    fwrite(fid, config.idx_neighbors_dense_grid, 'uint8');
    fwrite(fid, config.weights_neighbors_dense_grid, 'float');
    fwrite(fid, config.maxShift_dense, 'uint8');
    fwrite(fid, config.combination_ptr, 'int8');
    fwrite(fid, config.target_grid, 'float');                       % target grid
    
    fclose(fid);
end

