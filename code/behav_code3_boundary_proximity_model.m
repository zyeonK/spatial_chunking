
rng(918)
%% Set directory
data_dir = '../NC_data';
data_name  = 'sample_data.mat';
%% Load data
load(fullfile(data_dir,data_name)) % table_all

sbj_list = unique(table_all.sbj)';
num_of_sbj = length(sbj_list);

%% Boundary proximity model predicted response pattern

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Target data for boundary proximity model
target_list = [table_all.coin_x, table_all.coin_y];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% model parameter
arena_size = 1;
bnd_prox_c = arena_size/2;
softmax_tau = 0.1; 

step_size = 0.01;
x_coord_list = step_size:step_size:1-step_size;
y_coord_list = step_size:step_size:1-step_size;
tag_nan_x = [1:5, 95:99]; 
tag_nan_y = [1:5, 95:99];

% boundary proximity model
func_model = @(x,y) 1 ./ ([x, arena_size-x, y, arena_size-y] + bnd_prox_c); % boundary proximity model
func_softmax = @(mat,tau) exp(mat/tau) / sum(exp(mat(:)/tau));
pnt_list = combvec(x_coord_list, y_coord_list)';
vec_list = arrayfun(@(i) func_model(pnt_list(i,1), pnt_list(i,2)) , 1:length(pnt_list), 'uni', 0);

% code
sbj_mat_final = {}; sbj_mat_final_znorm = {};
for sbj_i = sbj_list
    sbj_i
    sbj_tag = table_all.sbj == sbj_i;

    curr_cue_pos = target_list(sbj_tag, :);
    cue_vec_list = arrayfun(@(i) func_model(curr_cue_pos(i,1), curr_cue_pos(i,2)) , 1:length(curr_cue_pos), 'uni', 0);

    sbj_mat_all = [];
    parfor cue_i = 1:length(cue_vec_list)
        vec_cue = cue_vec_list{cue_i};
    
        dist_to_pnt = cellfun(@(x) norm(vec_cue - x), vec_list);
        dist_to_pnt = max(dist_to_pnt) - dist_to_pnt;
        dist_to_pnt = func_softmax(dist_to_pnt, softmax_tau);
        dist_to_pnt = reshape(dist_to_pnt, [length(x_coord_list), length(y_coord_list)]);
        sbj_mat_all(:,:,cue_i) = dist_to_pnt;
    end

    avg_map_sbj = mean(sbj_mat_all,3);
    avg_map_sbj(tag_nan_x, :) = NaN;   
    avg_map_sbj(:,tag_nan_y) = NaN; 

    sbj_bnd_prox_resp{sbj_i} = avg_map_sbj;
    sbj_bnd_prox_resp_znorm{sbj_i} = (avg_map_sbj - nanmean(avg_map_sbj(:)))/std(avg_map_sbj(:),'omitnan');
end

%% Visualization - individual pattern (predicted by boundary proximity model)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data for visualization
visu_data = sbj_bnd_prox_resp_znorm;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

num_fig_col = 6;
num_fig_row = 3;

figure;
hold on
for err_i = 1:length(visu_data)
    subplot(num_fig_row, num_fig_col, err_i)
    hold on

    heatmap_data = visu_data{err_i}';
    imagesc(heatmap_data, 'AlphaData',~isnan(heatmap_data)); 
    set(gca,'YDir','normal')
    colormap('jet'); 
%     c = colorbar; 
%     c.LineWidth = 1.25;

    axis equal;
    clim([-3 3])
    set(gca,'LineWidth',1.5, 'FontWeight', 'bold','FontSize',15);
    title(['Sbj #', num2str(err_i)],'FontName','Helvetica','FontSize',17,'FontWeight','bold')
    box on
    % xticks([0+0.5, size(heatmap_data,1)/2+0.5, size(heatmap_data,1)+0.5])
    % xticklabels({"0", "0.5", "1"})
    % yticks([size(heatmap_data,1)/2+0.5, size(heatmap_data,1)+0.5])
    % yticklabels({"0.5", "1"})
    xticklabels({}); yticklabels({})    
    xlim([0.5 size(heatmap_data,1)+0.5])
    ylim([0.5 size(heatmap_data,1)+0.5])

end

%% Group-level response pattern (predicted by boundary proximity model)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data for group-level analysis
use_data = sbj_bnd_prox_resp_znorm; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% average pattern
avg_data = cat(3, use_data{:}); 
avg_data = mean(avg_data, 3,'omitnan'); 

% visualization - average pattern
figure;
hold on
heatmap_data = avg_data';
imagesc(heatmap_data, 'AlphaData',~isnan(heatmap_data)); %%%%%%%%
set(gca,'YDir','normal')
colormap('jet'); 
axis equal;
% c = colorbar; 
% c.LineWidth = 1.25;
clim([-1.5 1.5])
set(gca,'LineWidth',1.5, 'FontWeight', 'bold','FontSize',15);
box on
xticks([0+0.5, size(heatmap_data,1)/2+0.5, size(heatmap_data,1)+0.5])
xticklabels({"0", "0.5", "1"})
yticks([size(heatmap_data,1)/2+0.5, size(heatmap_data,1)+0.5])
yticklabels({"0.5", "1"})
xlim([0.5 size(heatmap_data,1)+0.5])
ylim([0.5 size(heatmap_data,1)+0.5])


% t-test for each pixel
tmap_all = [];
pval_all = [];
for row_i = 1:size(use_data{1},1)
    for col_i = 1:size(use_data{1},2)
        if sum(isnan(cellfun(@(x) x(row_i, col_i), use_data(:)))) == length(use_data)-1
            tmap_all(row_i, col_i) = NaN;
            pval_all(row_i, col_i) = NaN;
            continue;
        end
        [h,p,ci,stats] = ttest(cellfun(@(x) x(row_i, col_i), use_data(:)));
        tmap_all(row_i, col_i) = stats.tstat;
        pval_all(row_i, col_i) = p;
    end
end

% visualization - t-map (Fig 1d)
figure;
hold on
heatmap_data = tmap_all';
imagesc(heatmap_data, 'AlphaData',~isnan(heatmap_data)); 
set(gca,'YDir','normal')
colormap('jet'); 
c = colorbar; 
axis equal;
c.LineWidth = 1.25;
clim([-6 6])
set(gca,'LineWidth',1.2, 'FontWeight', 'bold','FontSize',15);
box on
xticks([0+0.5, size(heatmap_data,1)/2+0.5, size(heatmap_data,1)+0.5])
xticklabels({"0", "0.5", "1"})
yticks([size(heatmap_data,1)/2+0.5, size(heatmap_data,1)+0.5])
yticklabels({"0.5", "1"})
xlim([0.5 size(heatmap_data,1)+0.5])
ylim([0.5 size(heatmap_data,1)+0.5])


