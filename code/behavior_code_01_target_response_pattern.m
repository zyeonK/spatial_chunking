clear all
rng(918)

%% Set directory
data_dir = '../data';
data_name = 'sample_data.mat';

%% Load data
load(fullfile(data_dir,data_name)) % table_all
addpath('visualization_toolbox')

sbj_list = unique(table_all.sbj)';
num_of_sbj = length(sbj_list);

num_fig_col = 6;
num_fig_row = 3;

%% Target coin/ response (scatter plot)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data for visualization
visu_data = [table_all.resp_x, table_all.resp_y]; % response
% visu_data = [visu_table.coin_x,visu_table.coin_y] % target coin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Supplementary Fig 1a
figure;
hold on
scatter(visu_data(:,1), visu_data(:,2), 10, 'MarkerFaceColor','k', 'MarkerFaceAlpha',0.5, 'MarkerEdgeColor','none')
axis equal
box on
set(gca,'LineWidth',1.5, 'FontWeight', 'normal','FontSize',15);
xlim([0 1]); ylim([0 1])
xticks([0 0.5 1]); yticks([0.5 1])
xticklabels({}); yticklabels({})


% Supplementary Fig 1b
figure;
for sbj_i = 1:length(sbj_list) 
    tag_sbj = table_all.sbj == sbj_list(sbj_i);
    curr_data = visu_data(tag_sbj,:);

    subplot(num_fig_row, num_fig_col, sbj_i)
    hold on
    scatter(curr_data(:,1), curr_data(:,2), 10, 'MarkerFaceColor','k', 'MarkerFaceAlpha',0.5, 'MarkerEdgeColor','none')
    axis equal
    box on
    set(gca,'LineWidth',1.5, 'FontWeight', 'normal','FontSize',15);
    xlim([0 1]); ylim([0 1])
    xticks([0 0.5 1]); yticks([0.5 1])
    xticklabels({}); yticklabels({})
end

%% Target coin/ response pattern (heatmap)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data for generating pattern
use_data = table_all; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parameters
step_size = 0.01;
win_size = 0.1;
gaussianweight_sigma = 0.05;
x_coord_list = step_size:step_size:1-step_size;
y_coord_list = step_size:step_size:1-step_size;
tag_nan_x = [1:5, 95:99]; 
tag_nan_y = [1:5, 95:99];
            
% code
sbj_coin_all = {}; sbj_resp_all = {};
parfor sbj_i = 1:length(sbj_list) 

    sbj_tag = use_data.sbj == sbj_list(sbj_i);
    sbj_table = use_data(sbj_tag,:); 
    coin_x = sbj_table.coin_x;
    coin_y = sbj_table.coin_y;
    resp_x = sbj_table.resp_x;
    resp_y = sbj_table.resp_y;

    num_coin = []; num_resp = []; num_coin_gau = []; num_resp_gau = []; err_coin = []; score_coin = [];
    for x_coord_idx = 1:length(x_coord_list)
        x_coord = x_coord_list(x_coord_idx);
    
        for y_coord_idx = 1:length(y_coord_list)
            y_coord = y_coord_list(y_coord_idx);
    
            curr_pnt_coord = [x_coord, y_coord];

            coin_to_pnt = curr_pnt_coord - [coin_x, coin_y];
            resp_to_pnt = curr_pnt_coord - [resp_x, resp_y];

            % calculate distance between point and coin/response
            coin_to_pnt_dist = hypot(coin_to_pnt(:,1), coin_to_pnt(:,2));
            resp_to_pnt_dist = hypot(resp_to_pnt(:,1), resp_to_pnt(:,2));

            % gaussian weight based on the distance calculated above
            gaussianweight_coin = exp(-((coin_to_pnt_dist.^2)./(2*(gaussianweight_sigma.^2))));
            gaussianweight_resp = exp(-((resp_to_pnt_dist.^2)./(2*(gaussianweight_sigma.^2))));
            
            % tag within window coin/response
            tag_window_coin = ((x_coord-win_size/2)<=coin_x)&((x_coord+win_size/2)>=coin_x)&((y_coord-win_size/2)<=coin_y)&((y_coord+win_size/2)>=coin_y);
            tag_window_resp = ((x_coord-win_size/2)<=resp_x)&((x_coord+win_size/2)>=resp_x)&((y_coord-win_size/2)<=resp_y)&((y_coord+win_size/2)>=resp_y);

            % number of coin/response within current window
            num_coin(x_coord_idx, y_coord_idx) = sum(tag_window_coin);
            num_resp(x_coord_idx, y_coord_idx) = sum(tag_window_resp);


            % gaussian weighted
            num_coin_gau(x_coord_idx, y_coord_idx) = sum(tag_window_coin.*gaussianweight_coin);
            num_resp_gau(x_coord_idx, y_coord_idx) = sum(tag_window_resp.*gaussianweight_resp);

            if ~(sum(tag_window_resp)) && ~(sum(tag_window_coin)) % if there is no response or coin --> NaN
                num_resp(x_coord_idx, y_coord_idx) = NaN;
                num_resp_gau(x_coord_idx, y_coord_idx) = NaN;
            end
        end
    end


    num_coin(tag_nan_x, :) = NaN;
    num_coin(:,tag_nan_y) = NaN; 

    num_resp(tag_nan_x, :) = NaN;
    num_resp(:,tag_nan_y) = NaN; 

    num_resp_gau(tag_nan_x, :) = NaN;
    num_resp_gau(:,tag_nan_y) = NaN; 

    num_coin_gau(tag_nan_x, :) = NaN;
    num_coin_gau(:,tag_nan_y) = NaN; 

    % z-normalization within subject
    num_coin_gau = (num_coin_gau - nanmean(num_coin_gau(:)))./std(num_coin_gau(:),'omitnan');
    num_resp_gau = (num_resp_gau - nanmean(num_resp_gau(:)))./std(num_resp_gau(:),'omitnan');

    sbj_coin_all{sbj_i} = num_coin_gau; %num_coin;  
    sbj_resp_all{sbj_i} = num_resp_gau; %num_resp; 

end

% % save data
% save([data_dir, '/sbj_pattern_all.mat'],'sbj_coin_all','sbj_resp_all')

%% Visualization - individual pattern

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data for visualization
visu_data = sbj_resp_all; %sbj_coin_all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Supplementary Fig 1c
figure;
hold on
for err_i = 1:length(visu_data)
    subplot(num_fig_row, num_fig_col, err_i)
    hold on

    heatmap_data = visu_data{err_i}';
    imagesc(heatmap_data, 'AlphaData',~isnan(heatmap_data)); 
    set(gca,'YDir','normal')
    colormap('jet'); 
    axis equal;
    clim([-3 3])
    set(gca,'LineWidth',1.5, 'FontWeight', 'bold');
    box on
    xticklabels({}); yticklabels({})    
    xlim([0.5 size(heatmap_data,1)+0.5])
    ylim([0.5 size(heatmap_data,1)+0.5])

end

%% Group-level target coin/ response pattern

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data for group-level analysis
use_data = sbj_resp_all; % sbj_coin_all; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


avg_data = cat(3, use_data{:}); 
avg_data = mean(avg_data, 3,'omitnan'); 

% t-test for each pixel
tmap_all = [];
pval_all = [];
for row_i = 1:size(use_data{1},1)
    for col_i = 1:size(use_data{1},2)
        if sum(isnan(cellfun(@(x) x(row_i, col_i), use_data(:)))) >= (length(use_data)-1)
            tmap_all(row_i, col_i) = NaN;
            pval_all(row_i, col_i) = NaN;
            continue;
        end
        
        [h,p,ci,stats] = ttest(cellfun(@(x) x(row_i, col_i), use_data(:)));
        tmap_all(row_i, col_i) = stats.tstat;
        pval_all(row_i, col_i) = p;
    end
end

% visualization - t-map (Fig 2a)
figure;
hold on
heatmap_data = tmap_all';
imagesc(heatmap_data, 'AlphaData',~isnan(heatmap_data)); 
set(gca,'YDir','normal')
colormap('jet'); 
axis equal;
clim([-6 6])
set(gca,'LineWidth',1.2, 'FontWeight', 'normal','FontSize',15);
box on
xticks([0+0.5, size(heatmap_data,1)/2+0.5, size(heatmap_data,1)+0.5])
xticklabels({"0", "0.5", "1"})
xticklabels({})
yticks([size(heatmap_data,1)/2+0.5, size(heatmap_data,1)+0.5])
yticklabels({"0.5", "1"})
yticklabels({})
xlim([0.5 size(heatmap_data,1)+0.5])
ylim([0.5 size(heatmap_data,1)+0.5])

% save data
save([data_dir, '/group_pattern_data.mat'],'avg_data','tmap_all')

%% Leave-one-out pattern similarity analysis

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data for pattern similarity analysis
use_data = sbj_resp_all; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

loo_ps_r_all = []; 
loo_ps_fishz_all = [];
for sbj_i = 1:length(use_data)
    if isnan(use_data{sbj_i})
        loo_ps_r_all(sbj_i) = NaN; 
        loo_ps_fishz_all(sbj_i) = NaN; 
        continue;
    end
    leave_one_data = cat(3, use_data{~(sbj_i == [1:length(use_data)])}); 
    leave_one_avg_data = mean(leave_one_data, 3,'omitnan'); 
    [sbj_rho, sbj_p] = corr(use_data{sbj_i}(:), leave_one_avg_data(:),'rows','pairwise'); % using average map

    loo_ps_r_all(sbj_i) = sbj_rho; 
    loo_ps_fishz_all(sbj_i) = atanh(sbj_rho); 
end

% t-test for correlation coefficients
[h,p,ci,stats] = ttest(loo_ps_fishz_all); %% loo_ps_r_all 
d = stats.tstat / sqrt(stats.df+1);
fprintf('mean = %.3f, STD = %.3f\n', nanmean(loo_ps_fishz_all), std(loo_ps_fishz_all,'omitnan'))
fprintf('p = %.3f, t(%d) = %.3f, CI = [%.3f %.3f], Cohen''s d = %.3f\n',p, stats.df, stats.tstat, ci(1),ci(2), d)

