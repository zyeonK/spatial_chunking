rng(918)

%% Set directory
data_dir = '../NC_data';
data_name  = 'sample_data.mat';
%% Load data
load(fullfile(data_dir,data_name)) % table_all

sbj_list = unique(table_all.sbj)';
num_of_sbj = length(sbj_list);

%% Error vector data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data for error vector
err_vec_table = table_all(table_all.err<outlier_thres,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% parameters
win_size = 0.1;
step_size = 0.05;
x_coord_list = [step_size:step_size:1-step_size];
y_coord_list = [step_size:step_size:1-step_size];
outlier_thres = 0.2;

coin_x = err_vec_table.coin_x;
coin_y = err_vec_table.coin_y;
resp_x = err_vec_table.resp_x;
resp_y = err_vec_table.resp_y;

% code
norm_dir_vector = {};
for x_coord_idx = 1:length(x_coord_list)-1
    x_coord = x_coord_list(x_coord_idx);
    for y_coord_idx = 1:length(y_coord_list)-1
        y_coord = y_coord_list(y_coord_idx);

        tag_window_target = ((x_coord-win_size/2)<=coin_x)&((x_coord+win_size/2)>=coin_x)&((y_coord-win_size/2)<=coin_y)&((y_coord+win_size/2)>=coin_y);
        
        curr_coin_loc = [coin_x(tag_window_target),coin_y(tag_window_target)];
        curr_resp_loc = [resp_x(tag_window_target),resp_y(tag_window_target)];

        if size(curr_coin_loc,1) < 7
            norm_dir_vector{x_coord_idx, y_coord_idx} = NaN;
            continue;
        end

        norm_dir_vector{x_coord_idx, y_coord_idx} = mean((curr_resp_loc - curr_coin_loc)./sqrt(sum((curr_resp_loc - curr_coin_loc).^2,2)),1);
    end
end

x_dir = cellfun(@(x) x(1), norm_dir_vector);
y_dir = cellfun(@(x) x(end), norm_dir_vector);
x_dir(:,[1,end]) = NaN;
x_dir([1,end],:) = NaN;


%% Visualization - error vector map
coord_arena = combvec([1:size(norm_dir_vector,1)],[1:size(norm_dir_vector,2)]);
x_coord = coord_arena(1,:);
y_coord = coord_arena(2,:);

figure;
quiver(x_coord', y_coord', x_dir(:), y_dir(:),'k','linewidth',1.5,'MaxHeadSize',100,'AutoScaleFactor',1.5)
axis equal
box on
set(gca,'LineWidth',1.2, 'FontWeight', 'bold','FontSize',15);
title('Error vector','FontName','Helvetica','FontSize',20, 'FontWeight','bold')
xticks([0+0.5, size(norm_dir_vector,1)/2+0.5, size(norm_dir_vector,1)+0.5])
xticklabels({"0", "0.5", "1"})
yticks([size(norm_dir_vector,2)/2+0.5, size(norm_dir_vector,2)+0.5])
yticklabels({"0.5", "1"})
xlim([.5 size(norm_dir_vector,1)+0.5])
ylim([.5 size(norm_dir_vector,2)+0.5])
