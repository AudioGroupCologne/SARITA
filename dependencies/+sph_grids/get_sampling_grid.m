function [grid_data, num_nodes, N_sg, radius] = get_sampling_grid( ...
                                    grid_type, N, type, convention, plot)
% Wrapper to generate different spherical grid data. Per default grid is
% returned in radiant and colatitudes. Depending on grid type, default
% radii are defined.
%
% Parameter:
% ------------
%    grid_type:     available grids: 'gauss', 'lebedev', 'fliege',
%                                    'eigenmike32', 'hosma', 'zylia'
%                                    'horizontal', 'extremal', 'equiangular'
%    N:             desired order N of the sampling scheme
%    type:          string 'rad' or 'deg' {default: deg}
%    convention:    'col' : theta ranging from 0 to 180, as it is used by
%                           miro, sofia, ...
%                   'el' : theta ranging from 90 to 0 to -90, as it is used
%                          by SOFA convention, ...
%                   {default: 'col'}
%    plot:          plot sampling scheme {default: false}
%
% Returns:
% -----------
%   grid_data: M X [az, col, weights] in rad {or deg if specified}
%              with M: number of sampling positions
%
%   num_nodes: number of sampling positions
%
%   N_sg: order of the sampling scheme
%
%   radius: radius of the spherical microphone array in case its specified
%           by the sampling scheme
%
% Dependencies:
% -------------
% SOFIA toolbox - http://audiogroup.web.th-koeln.de/SOFiA_wiki/WELCOME.html
% getEigenmikeNodes - https://de.mathworks.com/matlabcentral/fileexchange/73721-geteigenmikenodes
%
% References:
% -------------
%  grid data for extremal, equiangular, and fliege-maier grids are taken
%  from subdeq toolbox - https://github.com/AudioGroupCologne/SUpDEq
%  extremal: https://web.maths.unsw.edu.au/~rsw/Sphere/Extremal/New/index.html
%
% (C) 01/2020 Tim Luebeck
%     latest update: 07.05.2021 Hannes Helmholz

    %% define some constants/presets...
    % Lebedev data from SOFIA toolbox
    lbdv_order = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 17, 20, 23, 26, 29, 32, 35, 38, 41, 44, 47, 50, 53, 56, 59, 62, 65];
    lbdv_num_nodes = [6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302, 350, 434, 590, 770, ...
                      974, 1202, 1454, 1730, 2030, 2354, 2702, 3074, 3470, 3890, 4334, 4802, 5294, 5810];

    % zylia cartesian data from "https://www.zylia.co/white-paper.html"
    zylia_grid_cart = [  0.0,    0.0,  49.0; ...
                        32.7,    0.1,  36.5; ...
                        -16.4,  28.3,  36.5; ...
                        -16.3, -28.3,  36.5; ...
                        6.3,   -45.8,  16.3; ...
                        36.6,  -28.2,  16.3; ...
                        36.5,   28.4,  16.3; ...
                        6.2,    45.8,  16.3; ...
                        -42.8,  17.4,  16.3; ...
                        -42.7, -17.6,  16.3; ...
                        -36.5, -28.4, -16.3; ...
                        -6.2,  -45.8, -16.3; ...
                        42.8,  -17.4, -16.3; ...
                        42.7,   17.6, -16.3; ...
                        -6.3,   45.8, -16.3; ...
                        -36.6,  28.2, -16.3; ...
                        -32.7,  -0.1, -36.5; ...
                        16.4,  -28.3, -36.5; ...
                        16.3,   28.3, -36.5];
    %%
    if nargin < 2 || isempty(N) || N < 0
       if any(strcmpi(grid_type, {'gauss', 'lebedev', 'fliege', 'extremal', 'equiangular'}))
            N = 1;
            disp('WARNING: sampling grid order N is not specified, use N=1');
       elseif strcmpi(grid_type, 'eigenmike32')
            N = 4;
            disp('WARNING: eigenmike32 is for order N=4');
       elseif strcmpi(grid_type, 'zylia')
            N = 3;
            disp('WARNING: zylia is for order N=3');
       elseif strcmpi(grid_type, 'hosma')
            N = 7;
            disp('WARNING: hosma 7n MK2 is for order N=7');
       elseif any(strcmp(grid_type, {'EMA', 'horizontal'}))
            disp('WARNING: No number of sampling points specified for pure horizontal grid.\n Calculate 360 equidistant points along the horizontal plane.\n');
            N = 360;
       end
    end
    if nargin < 3 || isempty(type)
        type = 'rad';
    end
    if nargin < 4 || isempty(convention)
        convention = 'col';
    end
    if nargin < 5 || isempty(plot)
        plot = 0;
    end

    %%
    if strcmpi(grid_type, 'gauss')
        az_nodes = 2 * (N+1);
        el_nodes = az_nodes/2;
        [grid_data, num_nodes, N_sg] = sofia_gauss(az_nodes, el_nodes, 0);
        radius = [];

    elseif strcmpi(grid_type, 'lebedev')
        lbdv_degree = lbdv_num_nodes(lbdv_order == N);
        if ~isempty(lbdv_degree)
            [grid_data, num_nodes, N_sg] = sofia_lebedev(lbdv_degree, 0);
        else
            idx = find((lbdv_order-N)>0, 1) - 1;
            N_next = lbdv_order(idx);
            lbdv_degree_next = lbdv_num_nodes(idx);
            fprintf('WARNING: Lebedev grid of order N=%d is not defined, falling back to N=%d\n', N, N_next);
            [grid_data, num_nodes, N_sg] = sofia_lebedev(lbdv_degree_next, 0);
        end
        radius = [];

    elseif strcmpi(grid_type, 'fliege')
        if N > 29
            disp('WARNING: Sorry, we just implemented Fliege grid orders up to 29.')
            N = 29;
        end
        load('fliegeMaierNodes_1_30.mat', 'fliegeNodes');
        grid_cart = fliegeNodes{N+1}(:, 1:3);
        [grid_data(:, 1), grid_data(:, 2)] = ...
                cart2sph(grid_cart(:, 1), grid_cart(:, 2), grid_cart(:, 3));
        % make azimuth ranging from zero to 2 pi
        grid_data(:, 1) = mod(grid_data(:, 1) + 2*pi, 2*pi);
        % make colatitudes from elevations
        grid_data(:, 2) = pi/2 - grid_data(:, 2);
        % calc weights
        grid_data(:, 3) = fliegeNodes{N+1}(:, 4) ./ sum(fliegeNodes{N+1}(:, 4));

        num_nodes = size(grid_data, 1);
        N_sg  = N;
        radius = [];

    elseif strcmpi(grid_type, 'extremal')
        if N > 100
            disp('WARNING: Sorry, so far we just implemented extremal grid orders up to 100.')
            N = 100;
        end
        load('extremalNodes_1_100.mat', 'extremalNodes');
        grid_data = deg2rad(extremalNodes{N}(:, 1:2));
        % get weights
        grid_data(:, 3) = extremalNodes{N}(:, 3);

        num_nodes = size(grid_data, 1);
        N_sg  = N;
        radius = [];

    elseif strcmpi(grid_type, 'equiangular')
        if N > 44
            disp('WARNING: Sorry, so far we just implemented equiangular grid orders up to 44.')
            N = 44;
        end
        load('equiangularNodes_1_44.mat', 'equiangularNodes');
        grid_data = deg2rad(equiangularNodes{N}(:, 1:2));
        % calc weights
        grid_data(:, 3) = 1/size(grid_data, 1) .* ones(size(grid_data, 1), 1);
        disp('WARNING: Eqiangular grids dont have quadrature weights.')

        num_nodes = size(grid_data, 1);
        N_sg  = N;
        radius = [];

    elseif strcmpi(grid_type, 'eigenmike32')
        if ~isempty(N) && N ~= 4
            disp('WARNING: Eigenmike just supports 32 sampling point grid, for SH processing up to order 4.')
        end
        grid_data = getEigenmikeNodes('rad', 0);
        grid_data(:, 2) = pi/2 - grid_data(:, 2);
        num_nodes = 32;
        N_sg = 4;
        radius = 0.042;

    elseif strcmpi(grid_type, 'hosma')
        if ~isempty(N) && N ~= 7
            disp('WARNING: Hosma just supports 64 sampling point grid for SH processing up to order 7.')
        end
        load('hosma_gridpoints.mat', 'fliege_sph_rot');
        grid_data = fliege_sph_rot(:, [1, 2, 4]); % kill radius
        grid_data(:, 1:2) = deg2rad(grid_data(:, 1:2));
        % make colatitudes from elevations
        grid_data(:, 2) = pi/2 - grid_data(:, 2);

        num_nodes = 64;
        N_sg = 7;
        radius = 0.11;

    elseif strcmpi(grid_type, 'zylia')
        if ~isempty(N) && N ~= 3
            disp('WARNING: Zylia just supports 19 sampling point grid for SH processing up to order 3.')
        end
        % matlab x front, y left, z top
        % zylia  x
        grid_data_cart(:, 1) = zylia_grid_cart(:, 1);
        grid_data_cart(:, 2) = zylia_grid_cart(:, 2);
        grid_data_cart(:, 3) = zylia_grid_cart(:, 3);

        [grid_data(:, 1), grid_data(:, 2), grid_data(:, 3)]  = ...
            cart2sph(grid_data_cart(:, 1), grid_data_cart(:, 2), grid_data_cart(:, 3) );
        grid_data(:, 1) = mod(grid_data(:, 1) + 2*pi, 2*pi);
        % make colatitudes from elevations
        grid_data(:, 2) = pi/2 - grid_data(:, 2);
        % weights
        grid_data(:, 3) = 1/size(grid_data, 1) .* ones(size(grid_data, 1), 1);
        disp('WARNING: Zylia grid dont has specified quadrature weights.')

        num_nodes = 19;
        N_sg = 3;
        radius = 0.049;
       
    elseif any(strcmpi(grid_type, {'EMA', 'horizontal'}))
        grid_data = [deg2rad([0 : 360/N : 360 - 360/N].'), ones(N, 1) * pi/2];
        num_nodes = size(grid_data, 1);
        N_sg = floor((num_nodes-1)/2);
        warning('Calculate sampling scheme order of an EMA')
    else
        error('Grid type not implemented so far.')
    end

    %%
    if strcmpi(convention, 'el')
        % Convert to elevation
        grid_data(:, 2) = pi/2 - grid_data(:, 2);
    end
    if strcmpi(type, 'deg')
        % Convert from rad 2 deg
        grid_data(:, 1:2) = grid_data(:, 1:2) * 180/pi;
    end
    if plot
       if strcmpi(type, 'deg')
           grid_data_2plot(:, 1:2) = grid_data(:, 1:2) * pi/180;
           if strcmpi(convention, 'el')
               grid_data_2plot(:, 2) = pi/2 - grid_data_2plot(:, 2);
           end
           plot_sampling_scheme(grid_data_2plot(:, 1:2))
       else
           if strcmpi(convention, 'el')
               grid_data_2plot(:, 1) = grid_data(:, 1);
               grid_data_2plot(:, 2) = pi/2 - grid_data(:, 2);
           else
               grid_data_2plot = grid_data;
           end
           plot_sampling_scheme(grid_data_2plot(:, 1:2))
       end
   end
end
