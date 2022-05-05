% This demo applies sarita upsampling on a simulated soundfield
% of single plane waves
%
%
%
% Tim Lübeck 03.11.2020
% latest update tim lübeck 04.01.2022

clear all; close all; clc

addpath('../../SARITA_dev/Sarita_code');%addpath(genpath('../Sarita_code'));
addpath(genpath('../dependencies'));
addpath(genpath('../src'))

do_plots = 1;
do_sound = 1;

disp('SARITA DEMO 1) Upsampling of a simulated Plane wave.')

%% (1) Define some parameters

array_config    = 2;      % rigid sphere
soft_limit  	= 20;     % radial filter softlimit
fs              = 48000;  % Sampling Frequency
NFFT        	= 512;    % number of FFT bins
c               = 343;
look_directions = pi/2; % for binaural signal calculation
% define sparsely measured grid
N_sparse = 5;
grid_data_sparse = sph_grids.get_sampling_grid('lebedev', N_sparse);
radius = 0.1;

% define reference/upsampling grid
N_ref = 29;
grid_data_ref = sph_grids.get_sampling_grid('lebedev', N_ref);

% estimate transit frequency 
f_transit = N_sparse * c/2/pi/radius*0.5;

%% (2) Simulate plane wave sound field of sparse and reference dataset

% direction of the first plane wave
az    = 0;           % azimuth angle
col   = pi/2;        % colatitude angle

time_delay    = 64/48000;    % avoids wrap arounds in binaural synthesis
delta_time    = 0.00082353;  % time shift between two plane waves
num_pws = 1;
atten = sqrt(0.5);

DRTFs_sparse  = zeros(size(grid_data_sparse, 1), round(NFFT/2) + 1);
DRTFs_ref     = zeros(size(grid_data_ref, 1), round(NFFT/2) + 1);

for pw_idx = 1:num_pws

    [DRTF_tmp, ~] = sofia_swg(radius, grid_data_sparse, array_config, fs, NFFT, az, col, 100, time_delay);
    DRTFs_sparse = DRTFs_sparse + DRTF_tmp .* atten^pw_idx;

    [DRTF_tmp, ~] = sofia_swg(radius, grid_data_ref, array_config, fs, NFFT, az, col, 100, time_delay);
    DRTFs_ref = DRTFs_ref + DRTF_tmp .* atten^pw_idx;

    % increase params for the next wave
    az = az + 30*(pi/180);
    col = col + 30*(pi/180);
    time_delay = time_delay + delta_time;

    clear DRTF_tmp
end

% inverse dft 
drirs = real(ifft([DRTFs_sparse, conj(fliplr(DRTFs_sparse(:, 2:end-1)))], [], 2)); 
drirs_ref = real(ifft([DRTFs_ref, conj(fliplr(DRTFs_ref(:, 2:end-1)))], [], 2)); 

%% (3) Perform upsampling of the sparse DRTF dataset
[drirs_ups, DRTFs_ups] = Sarita_upsampling(drirs, grid_data_sparse, grid_data_ref, radius, ...
                                             'f_transit', f_transit, 'frame_length', 32);


%% (4) Perform Spatial Fourier Transform
kr = (linspace(0, fs/2, NFFT/2 + 1) * 2 * pi * radius) / 343;  

DRTFs_sparse_nm = sofia_stc(N_sparse, DRTFs_sparse, grid_data_sparse);
DRTFs_ref_nm    = sofia_stc(N_ref, DRTFs_ref, grid_data_ref);
DRTFs_ups_nm    = sofia_stc(N_ref, DRTFs_ups, grid_data_ref);

radial_filters_sparse = sofia_mf(N_sparse, kr, array_config, soft_limit);
radial_filters_ref    = sofia_mf(N_ref, kr, array_config, soft_limit);

%% (5) Calculate binaural signals

[B_l, B_r] = sofia_binauralX(DRTFs_sparse_nm, radial_filters_sparse, look_directions);
brir_sparse = ifft(cat(3, [B_l, conj(B_l(:, end-1:-1:2, :))], ...
                          [B_r, conj(B_r(:, end-1:-1:2, :))]),  [], 2, 'symmetric');

[B_l, B_r] = sofia_binauralX(DRTFs_ref_nm, radial_filters_ref, look_directions);
brir_ref = ifft(cat(3, [B_l, conj(B_l(:, end-1:-1:2, :))], ...
                       [B_r, conj(B_r(:, end-1:-1:2, :))]),  [], 2, 'symmetric');
                 
[B_l, B_r] = sofia_binauralX(DRTFs_ups_nm, radial_filters_ref, look_directions);
brir_ups = ifft(cat(3, [B_l, conj(B_l(:, end-1:-1:2, :))], ...
                  [B_r, conj(B_r(:, end-1:-1:2, :))]),  [], 2, 'symmetric');

NFFT = max([size(brir_sparse, 2), size(brir_ref, 2), size(brir_ups, 2)]);

BRTFs_sparse = fft(brir_sparse, NFFT, 2); BRTFs_sparse = BRTFs_sparse(:, 1:NFFT/2 +1);
BRTFs_ref = fft(brir_ref, NFFT, 2); BRTFs_ref = BRTFs_ref(:, 1:NFFT/2 +1);
BRTFs_ups = fft(brir_ups, NFFT, 2); BRTFs_ups = BRTFs_ups(:, 1:NFFT/2 +1);

%% (6) PLOTS
if do_plots
    % plot binaural signals
    fig1 = figure(1);
    subplot(2, 3, 1)
        plot(squeeze(brir_sparse(1, :, 1)), 'Color', [139/255, 0, 0])
        hold on;
        plot(squeeze(brir_sparse(1, :, 2)), 'Color', [0, 139/255, 139/255])
        legend('left', 'right')
        title('Sparse')
    subplot(2, 3, 2)
        plot(squeeze(brir_ref(1, :, 1)), 'Color', [139/255, 0, 0])
        hold on;
        plot(squeeze(brir_ref(1, :, 2)), 'Color', [0, 139/255, 139/255])
        legend('left', 'right')
        title('Ref')
    subplot(2, 3, 3)
        plot(squeeze(brir_ups(1, :, 1)), 'Color', [139/255, 0, 0])
        hold on;
        plot(squeeze(brir_ups(1, :, 2)), 'Color', [0, 139/255, 139/255])
        legend('left', 'right')
        title('Upsampled')

    subplot(2, 3, [4, 5, 6])
        semilogx(linspace(5*eps, fs/2, size(BRTFs_sparse, 2)), 20*log10(abs(squeeze(BRTFs_sparse(1, :, 1)))), 'Color', [0.8, 0.8, 0.8], 'LineWidth', 3)
        hold on;
        semilogx(linspace(5*eps, fs/2, size(BRTFs_ref, 2)), 20*log10(abs(squeeze(BRTFs_ref(1, :, 1)))), 'Color', [139/255, 0, 0],'LineWidth', 2)
        semilogx(linspace(5*eps, fs/2, size(BRTFs_ups, 2)), 20*log10(abs(squeeze(BRTFs_ups(1, :, 1)))), 'Color', [0, 139/255, 139/255], 'LineWidth', 2)

        ylim([-60, 5])
        xlim([20 20000])
        grid on;
        legend('Sparse', 'Reference', 'Upsampled')
    
    %%
    % plot SOFIA MTX data
    f_observe = 15000;
    frequencies = linspace(5*eps,fs/2, size(DRTFs_sparse_nm, 2));
    
    [~ , krIndex] = min(abs(frequencies - f_observe));

    mtxData_sparse = sofia_makeMTX(N_ref, DRTFs_sparse_nm, radial_filters_ref, krIndex);
    mtxData_upsampled = sofia_makeMTX(N_ref, DRTFs_ups_nm, radial_filters_ref, krIndex);
    mtxData_ref = sofia_makeMTX(N_ref, DRTFs_ref_nm, radial_filters_ref, krIndex);

    Visualization_Style= 0; %Sphere
    %Visualization_Style= 1; %Flat

    fig2 = figure(2);
    clf();
    sofia_visual3D(mtxData_sparse, Visualization_Style);

    fig3 = figure(3);
    clf();
    sofia_visual3D(mtxData_upsampled, Visualization_Style);

    fig4 = figure(4);
    clf();
    sofia_visual3D(mtxData_ref, Visualization_Style);
    pause(0.01)


    % energy per order plots
    shift_value = 10*log10(N_ref/N_sparse);
    [~, ~, ~, avgEperN_Pnm] = AKshEnergy(DRTFs_sparse_nm);
    [~, ~, ~, avgEperN_Pnm_upsampled] = AKshEnergy(DRTFs_ups_nm);
    [~, ~, ~, avgEperN_Pnm_ref] = AKshEnergy(DRTFs_ref_nm);

    figure
    plot(0:N_ref,10*log10(avgEperN_Pnm_ref),'k','LineWidth',1.5)
    hold on
    plot(0:N_sparse,10*log10(avgEperN_Pnm)-shift_value,'r','LineWidth',1.5)
    plot(0:N_ref,10*log10(avgEperN_Pnm_upsampled),'b','LineWidth',1.5)
    legend('Ref','Sparse','Upsampled');
    title('Total Energy per Order N')
    xlabel('Order N')
    ylabel('Energy in dB')
    ylim([round(min(10*log10(avgEperN_Pnm_ref)))-2 round(max(10*log10(avgEperN_Pnm)-shift_value))+2]);
end

%% (8) - Optional: Listen to the results
if do_sound
    % load drum test signal
    [testSignal, ~] = audioread('MARA_LE_DRUMS.wav');
    clear bin_signal_*

    bin_signal_sparse(:, 1) = conv(testSignal(200:200+ 5*fs, :), brir_sparse(1, :, 1));
    bin_signal_sparse(:, 2) = conv(testSignal(200:200+ 5*fs, :), brir_sparse(1, :, 2));

    bin_signal_ref(:, 1) = conv(testSignal(200:200+ 5*fs, :), brir_ref(1, :, 1));
    bin_signal_ref(:, 2) = conv(testSignal(200:200+ 5*fs, :), brir_ref(1, :, 2));

    bin_signal_ups(:, 1) = conv(testSignal(200:200+ 5*fs, :), brir_ups(1, :, 1));
    bin_signal_ups(:, 2) = conv(testSignal(200:200+ 5*fs, :), brir_ups(1, :, 2));
    
    %%
    clear sound;
    soundsc(bin_signal_sparse, fs);
    %%
    clear sound;
    soundsc(bin_signal_ref, fs);
    %%
    clear sound;
    soundsc(bin_signal_ups, fs);
end
