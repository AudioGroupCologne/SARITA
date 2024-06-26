function [N_grid, weights, name] = find_grid_type(grid_data)
 which_grids_to_compare = {'lebedev', 'gauss', 'fliege', 'extremal', 'equiangular', 'eigenmike32', 'hosma', 'zylia'};
 which_N_to_compare = [1:44];
 
 thresh = 1.0e-6;
 found = 0;
 
 name_idx = 1;
 
 N_grid = [];
 weights = [];
 name = [];

 while name_idx <= length(which_grids_to_compare) && ~found
     n_idx = 1;
     while n_idx <= length(which_N_to_compare) && ~found
         N = which_N_to_compare(n_idx);
         [grid_data_comp, ~, N_sg]  = sph_grids.get_sampling_grid(string(which_grids_to_compare{name_idx}), N);
       
         if size(grid_data_comp, 1) == size(grid_data, 1)
            check = abs(grid_data_comp(:, 1:2)-grid_data );
            if max(max(check)) < thresh
                N_grid = N_sg;
                weights = grid_data_comp(:, 3);
                name = string(which_grids_to_compare{name_idx});
                fprintf('---------------------------------------\n')
                fprintf('Found N=%d %s grid\n', N_sg, name);
                found = 1;
            end
         end
         n_idx = n_idx + 1;
     end
     name_idx = name_idx + 1;
 end
 
for name_idx = 1:length(which_grids_to_compare)
   for N = which_N_to_compare
       
   end
end

end