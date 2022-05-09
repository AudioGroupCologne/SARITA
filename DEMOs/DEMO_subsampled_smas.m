% Simulate sparse array measurements by spatial subsampling, apply Sarita
% upsampling, and compare binaural signals
%
% Tim Lübeck 03.11.2020
%
clear all; close all; clc
do_plots = 1;
do_sound = 1;
addpath('../../SARITA_dev/Sarita_code'); % addpath('../SARITA_code');
addpath(genpath('../dependencies'));
addpath(genpath('../src'))

disp('SARITA DEMO 2) Upsampling of measured SMA signals.')
%% (1) load Data

N_ref = 29;
[drirs, grid_data_ref, fs, radius, ~] = sfe_import_sofa('DRIR_CR1_VSA_1202RS_R.sofa', 1);
[brirs_dh, grid_data_brir, ~, ~, ~] = sfe_import_sofa('BRIR_CR1_KU_MICS_R.sofa');

%% (2) define some parameters 
array_config = 2;      % rigid sphere
soft_limit   = 20;     % radial filter softlimit
NFFT         = max([size(drirs, 2), size(brirs_dh, 2)]);    % number of FFT bins
c            = 343;
look_directions = deg2rad([[0:359]', ones(360, 1)*90]); %[pi/2, pi/2]; % for binaural signal calculation

% define sparsely measured grid
N_sparse = 3;
grid_data_sparse = sph_grids.get_sampling_grid('lebedev', N_sparse);

% estimate transit frequency with kr < N
f_transit = N_sparse * c/2/pi/radius * 0.5;  % if set, low-frequencies are      
                                             % 'linear' weighted (not separately for mag and phase)

%% (3) simulate sparse array measurement by subsampling
DRTFs_ref = fft(drirs, [], 2); DRTFs_ref = DRTFs_ref(:, 1:round(length(DRTFs_ref)/2 + 1));
DRTFs_ref_nm = sofia_stc(N_ref, DRTFs_ref, grid_data_ref);
DRTFs_sparse = sofia_itc(DRTFs_ref_nm, grid_data_sparse);
drirs_sparse = real(ifft([DRTFs_sparse, conj(fliplr(DRTFs_sparse(:, 2:end-1)))], [], 2)); 

%% (4) Perform upsampling of the sparse DRTF dataset
[~, DRTFs_ups] = Sarita_upsampling(drirs_sparse, grid_data_sparse, grid_data_ref(:, 1:2), radius, ...
                                   'f_transit', f_transit, 'frame_length', 32);
                                         
%% (5) PerformsSpatial Fourier transform
kr = (linspace(0, fs/2, NFFT/2 + 1) * 2 * pi * radius) / 343;  

DRTFs_sparse_nm = sofia_stc(N_sparse, DRTFs_sparse, grid_data_sparse);
DRTFs_ups_nm    = sofia_stc(N_ref, DRTFs_ups, grid_data_ref);

%% (6) Calculate binaural signals

radial_filters_sparse = sofia_mf(N_sparse, kr, array_config, soft_limit);
radial_filters_ref  = sofia_mf(N_ref, kr, array_config, soft_limit);

[BRTF_l, BRTF_r] = sofia_binauralX(DRTFs_sparse_nm, radial_filters_sparse, look_directions);
brirs_sparse = ifft(cat(3, [BRTF_l, conj(BRTF_l(:, end-1:-1:2, :))], ...
                           [BRTF_r, conj(BRTF_r(:, end-1:-1:2, :))]),  [], 2, 'symmetric');
                   
[BRTF_l, BRTF_r] = sofia_binauralX(DRTFs_ref_nm, radial_filters_ref, look_directions);
brirs_ref = ifft(cat(3, [BRTF_l, conj(BRTF_l(:, end-1:-1:2, :))], ...
                        [BRTF_r, conj(BRTF_r(:, end-1:-1:2, :))]),  [], 2, 'symmetric');
                       

[BRTF_l, BRTF_r] = sofia_binauralX(DRTFs_ups_nm, radial_filters_ref, look_directions);
brirs_ups = ifft(cat(3, [BRTF_l, conj(BRTF_l(:, end-1:-1:2, :))], ...
                        [BRTF_r, conj(BRTF_r(:, end-1:-1:2, :))]),  [], 2, 'symmetric');

clear BRTF_l BRTF_r
%% (7) PLOTS

if do_plots
    % plot binaural signals
    NFFT = 2^nextpow2(max([length(brirs_sparse), length(brirs_ref), length(brirs_ups)]));
    
    BRTFs_sparse = fft(brirs_sparse, NFFT, 2); BRTFs_sparse = BRTFs_sparse(:, 1:end/2+1, :);
    BRTFs_ref = fft(brirs_ref, NFFT, 2); BRTFs_ref = BRTFs_ref(:, 1:end/2+1, :);
    BRTFs_ups = fft(brirs_ups, NFFT, 2); BRTFs_ups = BRTFs_ups(:, 1:end/2+1, :);
    
    % normalize signals to 200 Hz bin
    bin = find(linspace(0, fs/2, NFFT) < 200, 1, 'last');
    BRTFs_ref = BRTFs_ref ./ mean(abs(BRTFs_ref(1, bin, :)), 3);
    BRTFs_sparse = BRTFs_sparse ./ mean(abs(BRTFs_sparse(1, bin, :)), 3);
    BRTFs_ups = BRTFs_ups ./ mean(abs(BRTFs_ups(1, bin, :)), 3);
    
    brirs_ref = cat(3, ifft([BRTFs_ref(:, :, 1), conj(BRTFs_ref(:, end-1:-1:2, 1))], [], 2, 'symmetric'), ...
                       ifft([BRTFs_ref(:, :, 2), conj(BRTFs_ref(:, end-1:-1:2, 2))], [], 2, 'symmetric'));
    brirs_sparse = ifft([BRTFs_sparse, conj(BRTFs_sparse(:, end-1:-1:2, :))], [], 2, 'symmetric');
    brirs_ups = ifft([BRTFs_ups, conj(BRTFs_ups(:, end-1:-1:2, :))], [], 2, 'symmetric');
    
    fig1 = figure(1);
    subplot(2, 3, 1)
        plot(squeeze(brirs_sparse(1, :, 1)), 'Color', [139/255, 0, 0])
        hold on;
        plot(squeeze(brirs_sparse(1, :, 2)), 'Color', [0, 139/255, 139/255])
        legend('left', 'right')
        title('Sparse')
    subplot(2, 3, 2)
        plot(squeeze(brirs_ref(1, :, 1)), 'Color', [139/255, 0, 0])
        hold on;
        plot(squeeze(brirs_ref(1, :, 2)), 'Color', [0, 139/255, 139/255])
        legend('left', 'right')
        title('Ref')
    subplot(2, 3, 3)
        plot(squeeze(brirs_ups(1, :, 1)), 'Color', [139/255, 0, 0])
        hold on;
        plot(squeeze(brirs_ups(1, :, 2)), 'Color', [0, 139/255, 139/255])
        legend('left', 'right')
        title('Upsampled')
        
    subplot(2, 3, [4, 5, 6])  
        
        % 3rd octave smoothing 
        BRTFs_sparse_smoothed = AKfractOctSmooth(squeeze(BRTFs_sparse(1, :, 1)).', 'amp', fs).';
        BRTFs_ref_smoothed = AKfractOctSmooth(squeeze(BRTFs_ref(1, :, 1)).', 'amp', fs).';
        BRTFs_ups_smoothed = AKfractOctSmooth(squeeze(BRTFs_ups(1, :, 1)).', 'amp', fs).';
        
        semilogx(linspace(5*eps, fs/2, size(BRTFs_sparse, 2)), 20*log10(abs(BRTFs_sparse_smoothed(1, :))), 'Color', [0.8, 0.8, 0.8], 'LineWidth', 3)
        hold on;
        semilogx(linspace(5*eps, fs/2, size(BRTFs_ref, 2)), 20*log10(abs(BRTFs_ref_smoothed(1, :))), 'Color', [139/255, 0, 0],'LineWidth', 2)
        semilogx(linspace(5*eps, fs/2, size(BRTFs_ups, 2)), 20*log10(abs(BRTFs_ups_smoothed(1, :))), 'Color', [0, 139/255, 139/255], 'LineWidth', 2)
        
        legend('Sparse', 'Reference', 'Upsampled')
        
        xlabel('Frequency in Hz')
        xticks([100 1000, 10000, 16000, 20000])
        xlim([20 20000])
        grid on;
end

%% (9) Listen to the results

if do_sound
    
    fs = 48000;
    % load drum test signal
    [testSignal, ~] = audioread('../src/test_signals/MARA_LE_DRUMS.wav');
    clear bin_signal_*
    
    bin_signal_sparse(:, 1) = conv(testSignal(200:200+ 5*fs, :), brirs_sparse(1, :, 1));
    bin_signal_sparse(:, 2) = conv(testSignal(200:200+ 5*fs, :), brirs_sparse(1, :, 2));
    
    bin_signal_ref(:, 1) = conv(testSignal(200:200+ 5*fs, :), brirs_ref(1, :, 1));
    bin_signal_ref(:, 2) = conv(testSignal(200:200+ 5*fs, :), brirs_ref(1, :, 2));
    
    bin_signal_ups(:, 1) = conv(testSignal(200:200+ 5*fs, :), brirs_ups(1, :, 1));
    bin_signal_ups(:, 2) = conv(testSignal(200:200+ 5*fs, :), brirs_ups(1, :, 2));
    
    bin_signal_dh(:, 1) = conv(testSignal(200:200+ 5*fs, :), brirs_dh(1, :, 1));
    bin_signal_dh(:, 2) = conv(testSignal(200:200+ 5*fs, :), brirs_dh(1, :, 2));
    
    %% 
    clear sound;
    soundsc(bin_signal_sparse, fs);
    %%
    clear sound;
    soundsc(bin_signal_ref, fs);
    %%
    clear sound;
    soundsc(bin_signal_ups, fs);
    %%
    clear sound;
    soundsc(bin_signal_dh, fs);
end

 

