clear all
rng(918)

%% Set directory
data_dir = '../data';
data_name  = 'sample_data.mat';
%% Load data
load(fullfile(data_dir,data_name)) % table_all

sbj_list = unique(table_all.sbj)';
num_of_sbj = length(sbj_list);

%% Error vector data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data for error vector
err_vec_table = table_all(table_all.score > 0, :);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parameters
win_size = 0.1;
step_size = 0.05;
x_coord_list = step_size:step_size:1-step_size;
y_coord_list = step_size:step_size:1-step_size;

coin_x = err_vec_table.coin_x;
coin_y = err_vec_table.coin_y;
resp_x = err_vec_table.resp_x;
resp_y = err_vec_table.resp_y;

map_u = nan(length(x_coord_list), length(y_coord_list)); 
map_v = nan(length(x_coord_list), length(y_coord_list)); 


for x_coord_idx = 1:length(x_coord_list)
    for y_coord_idx = 1:length(y_coord_list)
        
        x_coord = x_coord_list(x_coord_idx);
        y_coord = y_coord_list(y_coord_idx);
        
        tag_window_target = ((x_coord-win_size/2)<=coin_x)&((x_coord+win_size/2)>=coin_x)&((y_coord-win_size/2)<=coin_y)&((y_coord+win_size/2)>=coin_y);
        
        curr_coin = [coin_x(tag_window_target),coin_y(tag_window_target)];
        curr_resp = [resp_x(tag_window_target),resp_y(tag_window_target)];

        err_vecs = curr_resp - curr_coin;
        unit_vecs = err_vecs ./ hypot(err_vecs(:,1), err_vecs(:,2));
        mean_vec = mean(unit_vecs, 1);
        map_u(x_coord_idx, y_coord_idx) = mean_vec(1);
        map_v(x_coord_idx, y_coord_idx) = mean_vec(2);
    end
end
% near boundaries - NaNs
map_u([1, end], :) = NaN; map_u(:, [1, end]) = NaN;
map_v([1, end], :) = NaN; map_v(:, [1, end]) = NaN;

%% Visualization - error vector map (Fig 2d)

[X_grid, Y_grid] = meshgrid(x_coord_list, y_coord_list);
plot_u = map_u'; 
plot_v = map_v';

figure;
hold on;
quiver(X_grid, Y_grid, plot_u, plot_v, 'k', 'linewidth', 2, 'AutoScaleFactor', 1.8);
axis equal;
box on;
set(gca, 'LineWidth', 1.2, 'FontWeight', 'bold', 'FontSize', 15);
xlim([0 1]); 
ylim([0 1]);
xticks([0, 0.5, 1]); 
yticks([0, 0.5, 1]);

