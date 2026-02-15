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

%% Target coin/ response pattern across errors (heatmap)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data for visualization
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
            
thres_list = num2cell(0:25:100);
err_thres = cellfun(@(x) prctile(use_data.err, x), thres_list);

% code
err_coin_map = {}; err_resp_map = {};
for sbj_i = 1:length(sbj_list)
    curr_data_sbj = {};
    for thres_i = 1:length(thres_list)-1

        target_err_thres = err_thres; 
        err_low_bnd = target_err_thres(thres_i);
        err_up_bnd = target_err_thres(thres_i+1);

        curr_table = use_data((use_data.sbj == sbj_list(sbj_i))&(use_data.err >= err_low_bnd)&(use_data.err <= err_up_bnd),:); 
        coin_x = curr_table.coin_x;
        coin_y = curr_table.coin_y;
        resp_x = curr_table.resp_x;
        resp_y = curr_table.resp_y;
    
        num_coin = []; num_resp = []; num_coin_gau = []; num_resp_gau = []; 
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
    
        % z-normalization
        num_coin_gau = (num_coin_gau - nanmean(num_coin_gau(:)))./std(num_coin_gau(:),'omitnan');
        num_resp_gau = (num_resp_gau - nanmean(num_resp_gau(:)))./std(num_resp_gau(:),'omitnan');
        curr_coin_sbj{thres_i} = num_coin_gau;
        curr_resp_sbj{thres_i} = num_resp_gau;

    end
    err_coin_map{sbj_i} = curr_coin_sbj;
    err_resp_map{sbj_i} = curr_resp_sbj;
end

%% Group-level target coin/ response pattern

% Fig 2b
figure;
hold on
num_of_errs = length(thres_list)-1;
num_fig_col = num_of_errs;
num_fig_row = 1;

err_avg_resp_all = {};
for err_i = 1:num_of_errs
    subplot(num_fig_row, num_fig_col, err_i)
    hold on
    get_data = cellfun(@(x) x(err_i), err_resp_map);
    avg_map_data = cat(3, get_data{:}); 
    avg_map_data = mean(avg_map_data, 3,'omitnan'); 
    err_avg_resp_all{err_i} = avg_map_data;
    
    % t-test for each pixel
    tmap_all = [];
    pval_all = [];
    for row_i = 1:size(avg_map_data,1)
        for col_i = 1:size(avg_map_data,2)
            if sum(isnan(cellfun(@(x) x(row_i, col_i), get_data(:)))) >= length(get_data)-1
                tmap_all(row_i, col_i) = NaN;
                pval_all(row_i, col_i) = NaN;
                continue;
            end

            [h,p,ci,stats] = ttest(cellfun(@(x) x(row_i, col_i), get_data(:)));
            tmap_all(row_i, col_i) = stats.tstat;
            pval_all(row_i, col_i) = p;
        end
    end

    heatmap_data = tmap_all';
    imagesc(heatmap_data, 'AlphaData', ~isnan(heatmap_data));
    set(gca,'YDir','normal')
    colormap('jet'); 
    axis equal;
    clim([-4 4])
    set(gca,'LineWidth',1.2, 'FontWeight', 'normal','FontSize',14);
    box on
    xticks([0.5, size(heatmap_data,1)/2+0.5, size(heatmap_data,1)+0.5])
    xticklabels({"0", "0.5", "1"})
    xticklabels({})
    yticks([size(heatmap_data,1)/2+0.5, size(heatmap_data,1)+0.5])
    yticklabels({"0.5", "1"})
    yticklabels({})
    xlim([0.5 size(heatmap_data,1)+0.5])
    ylim([0.5 size(heatmap_data,1)+0.5])

end

%% Leave-one-out pattern similarity analysis (across errors)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data for pattern similarity analysis
use_data = err_resp_map; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Pattern similarity
loo_ps_r_all = []; loo_ps_fishz_all = [];
for sbj_i = 1:num_of_sbj
    get_data = use_data{sbj_i}; 
        loo_ps_r = []; loo_ps_fishz = [];
    for err_i = 1:length(get_data)
        avg_map_data = cat(3, get_data{~(err_i == [1:length(get_data)])}); 
        avg_map_data = mean(avg_map_data, 3,'omitnan'); 
        [loo_ps_r(err_i),p] = corr(get_data{err_i}(:), avg_map_data(:),'rows','complete');
        loo_ps_fishz(err_i) = atanh(loo_ps_r(err_i)); 
    end
    loo_ps_r_all = [loo_ps_r_all; loo_ps_r];
    loo_ps_fishz_all = [loo_ps_fishz_all; loo_ps_fishz];
end

% t-test for correlation coefficients (Table 2)
for err_i = 1:size(loo_ps_fishz_all,2)
    curr_err_data = loo_ps_fishz_all(:,err_i);
    [h,p,ci,stats] = ttest(curr_err_data); 
    d = stats.tstat / sqrt(stats.df+1);
    fprintf('mean = %.3f, STD = %.3f\n', nanmean(curr_err_data), std(curr_err_data,'omitnan'))
    fprintf('p = %.3f, t(%d) = %.3f, CI = [%.3f %.3f], Cohen''s d = %.3f\n\n',p, stats.df, stats.tstat, ci(1),ci(2), d)
end
