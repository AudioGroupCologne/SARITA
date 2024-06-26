function plot_sampling_scheme(grid_sph_rad)
% Parameter:
% ------------
% grid_sph_rad:  spherical coordinates with azimuth and colatitude!
% Returns:
% ------------
%
% Dependencies:
% -------------
%
% (C) 01/2020 Tim Luebeck

radius = 1.01;
grid_cart = zeros(size(grid_sph_rad, 1), 3);
[grid_cart(:, 1), grid_cart(:, 2), grid_cart(:, 3)] = ...
    sph2cart(grid_sph_rad(:, 1), pi/2 - grid_sph_rad(:, 2), radius);

figure();
colormap Gray;
if size(grid_cart(:, 1), 1)>1500
    plot3(grid_cart(:, 1), grid_cart(:, 2), grid_cart(:, 3),'marker','.','markerfacecolor','g','color','g','linestyle','none')
else
    for (node_idx = 1 : length(grid_cart))
        %plot3(grid_cart(node_idx, 1), grid_cart(node_idx, 2), grid_cart(node_idx, 3), ...
        %    'marker','o', 'markerfacecolor', 'g', 'color','g','linestyle','none'); hold on;
        plot3(grid_cart(node_idx, 1), grid_cart(node_idx, 2), grid_cart(node_idx, 3)); hold on;
        text(grid_cart(node_idx, 1), grid_cart(node_idx, 2), grid_cart(node_idx, 3), num2str(node_idx), 'Color', 'g');
    end
end
%quiver3(0, 0, 0, radius,      0,      0, 'Color', 'r', 'LineWidth', 1.5, 'MaxHeadSize', 1);
%quiver3(0, 0, 0,      0, radius,      0, 'Color', 'g', 'LineWidth', 1.5, 'MaxHeadSize', 1);
%quiver3(0, 0, 0,      0,      0, radius, 'Color', 'b', 'LineWidth', 1.5, 'MaxHeadSize', 1);
% %plot3(1, 1,ref(3),'x','MarkerSize',5);
line([0 radius+0.3],[0 0],[0 0], 'LineWidth', 2.1, 'Color', 'black');
text(radius+0.35, 0, 0, 'X', 'FontSize', 12);

line([0 0],[0 radius+0.3],[0 0], 'LineWidth', 2.1, 'Color', 'black');
text(0, radius+0.35, 0, 'Y', 'FontSize', 12);

line([0 0],[0 0],[0, radius+0.3], 'LineWidth', 2.1, 'Color', 'black');
text(0, 0, radius+0.35, 'Z', 'FontSize', 12);

axis off;
hold on;
grid off;
sphere;
axis equal;
rotate3d on;
light;
alpha(.8);
lighting phong;
camzoom(1.4);
hold off;
end