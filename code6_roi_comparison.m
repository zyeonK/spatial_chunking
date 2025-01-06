
%% collect and save 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
reg_setting_name = 'task_period_basic_add_keypress';
% reg_setting_name = 'node_coin_binary';
reg_setting_name = 'node_resp_coin_binary';
% reg_setting_name = 'node_coin_binary_exclude_task';

ctx = '';
% ctx = '_ctx';
% ctx = '_hpc_ap';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load(sprintf('DATA_ORGANIZED%s_%s.mat', ctx, reg_setting_name))



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
neural_metric = cellfun(@(x) cellfun(@(y) y, x,'uni',0), metric_all, 'uni', 0);
roi_name_list
con_name_list
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

COLOR_EM = [255/256,143/256,120/256];
COLOR_SPATIAL = [159/256,204/256,78/256];
COLOR_WORD = [122/256,164/256,255/256];


%% p table

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
size_mult_x = 80;
size_mult_y = 20;

comparison_type = 1; % overall
% comparison_type = 2; % across group

p_cutoff = 0.01;

sbj_flag = true(1,length(group));
sbj_flag = group==1;
% sbj_flag = group==0;

% sbj_flag = cov >= median(cov);
% sbj_flag = cov < median(cov);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data = [];
for roi_i = 1:length(roi_name_list)
    for con_i = 1:length(con_name_list)
        if comparison_type == 1
            temp = neural_metric{con_i}{roi_i}(sbj_flag);
            [~,p] = ttest(temp);
            temp = nanmean(temp);
        elseif comparison_type == 2
            temp1 = neural_metric{con_i}{roi_i}(group==1);
            temp0 = neural_metric{con_i}{roi_i}(group==0);
            [~,p] = ttest2(temp1,temp0);
            p = ranksum(temp1,temp0);
            temp = nanmean(temp1) - nanmean(temp0);
        end
        
        if temp<0
            p = -p;
        end
        data(roi_i,con_i) = p;
    end
end
figure('Position',[0,0, size(data,2)*size_mult_x+100, size(data,1)*size_mult_y+300]);
jh_p_table(data, roi_name_list, cellfun(@(x) replace(x,'_',' '),con_name_list, 'uni', 0), 0.05)
figure('Position',[700,0, size(data,2)*size_mult_x+100, size(data,1)*size_mult_y+300]);
jh_p_table(data, roi_name_list, cellfun(@(x) replace(x,'_',' '),con_name_list, 'uni', 0), 0.1)
% caxis([0 0.1])

%

%% set behavior metrics - OLD VERSION

%%%%%%%%%%%%%%%%%%%%%%%%%% refer ../code_em1_behavior/code5_metrics_correlation

% % accuracy and training metric
% aux_func_get_residual = @(mdl) mdl.Residuals.Raw;
% func_get_residual = @(x,y) aux_func_get_residual(fitlm(x(:),y(:)));
% 
% sess1_behav_spatial_list = {spatial_metric_sess_wise{1}, spatial_metric_sess_wise_25{1}, spatial_metric_sess_wise_50{1}, spatial_metric_sess_wise_75{1} };
% sess8_behav_spatial_list = {spatial_metric_sess_wise{8}, spatial_metric_sess_wise_25{8}, spatial_metric_sess_wise_50{8}, spatial_metric_sess_wise_75{8} };
% sess8_behav_spatial_raw_list = {spatial_metric_sess_wise_raw{8}, spatial_metric_sess_wise_raw{8}, spatial_metric_sess_wise_raw{8}, spatial_metric_sess_wise_raw{8} };
% training_behav_spatial_list = {spatial_metric_training, spatial_metric_training_25, spatial_metric_training_50, spatial_metric_training_75};
% 
% training_cov_behav_spatial_list = {};
% template = nan(1,length(sess1_behav_spatial_list{1}));
% for i = 1:length(training_behav_spatial_list)
%     training_cov_behav_spatial_list{i} = template;
%     training_cov_behav_spatial_list{i}(group==1) = ...
%         func_get_residual(sess1_behav_spatial_list{i}(group==1), training_behav_spatial_list{i}(group==1));
% end
% 
% sess8_behav_spatial_raw_name_list = {'sess8_spatial_raw','sess8_spatial_raw_25','sess8_spatial_raw_50','sess8_spatial_raw_75'};
% sess8_behav_spatial_name_list = {'sess8_spatial_acc','sess8_spatial_acc_25','sess8_spatial_acc_50','sess8_spatial_acc_75'};
% training_behav_spatial_name_list = {'training_spatial','training_spatial_acc_25','training_spatial_acc_50','training_spatial_acc_75'};
% training_cov_behav_spatial_name_list = {'training_cov_spatial','training_cov_spatial_acc_25','training_cov_spatial_acc_50','training_cov_spatial_acc_75'};
% 
% 
% % chunking index
% load('chunking_index.mat')
% chunking_idx1 = nan(1,length(sbj_list));
% chunking_idx2 = chunking_idx1;
% chunking_idx1(group~=0) = chunking_all(1,:);
% chunking_idx2(group~=0) = chunking_all(2,:);
% 
% chunking_idx1_cov_raw = func_get_residual(spatial_metric_sess_wise_raw{8}, chunking_idx1);
% chunking_idx2_cov_raw = func_get_residual(spatial_metric_sess_wise_raw{8}, chunking_idx2);
% chunking_idx1_cov_acc = func_get_residual(spatial_metric_sess_wise{8}, chunking_idx1);
% chunking_idx2_cov_acc = func_get_residual(spatial_metric_sess_wise{8}, chunking_idx2);
% 
% chunking_idx_list = {chunking_idx1, chunking_idx2, chunking_idx1_cov_raw, chunking_idx2_cov_raw, ...
%                      chunking_idx1_cov_acc, chunking_idx2_cov_acc };
% chunking_index_name_list = {'chunking_idx_trunc_pixcorr','chunking_idx_trunc_snr', ...
%                             'chunking_idx_trunc_pixcorr_cov_raw','chunking_idx_trunc_snr_cov_raw', ...
%                             'chunking_idx_trunc_pixcorr_cov_acc','chunking_idx_trunc_snr_cov_acc'};
% 
% % combine
% behav_metric_list = [sess8_behav_spatial_list, sess8_behav_spatial_raw_list, ...
%                     training_behav_spatial_list, training_cov_behav_spatial_list, ...
%                     chunking_idx_list ];
% behav_metric_name_list = [ sess8_behav_spatial_name_list, sess8_behav_spatial_raw_name_list, ...
%                            training_behav_spatial_name_list, training_cov_behav_spatial_name_list, ...
%                            chunking_index_name_list ];

% behav_metric_name_list = cellfun(@(x) ['resp_',x], behav_metric_name_list, 'uni', 0);  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% set behavior metrics - OLD VERSION

%%%%%%%%%%%%%%%%%%%%%%%%%% refer ../code_em1_behavior/code5_metrics_correlation

% chunking index
load('chunking_idx.mat')


behav_metric_list = arrayfun(@(i) data(i,:), 1:size(data,1), uni=0);
behav_metric_name_list = arrayfun(@(i) sprintf('chunking%d',i),1:size(data,1), uni=0);


%%%%%%%%%
aux_func_get_residual = @(mdl) mdl.Residuals.Raw;
func_get_residual = @(x,y) aux_func_get_residual(fitlm(x(:),y(:)));

cov = [0.2331    0.1889    0.2439    0.1441    0.2554    0.2228    0.1702    0.1293    0.2229    0.2708    0.3125    0.1944    0.2448    0.1179    0.1660    0.1191    0.2660    0.2268    0.1530    0.3579    0.1642    0.2118    0.1294    0.2207    0.1929    0.3283    0.0820    0.0791    0.2131    0.1945    0.1107    0.1771    0.1615];
cov = cov(group==1);

% cov = spatial_metric_training(group==1);
% behav_metric_list = cellfun(@(x) func_get_residual(cov,x), behav_metric_list, uni=0);

%% correlation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
con_i = 1;


flag = group==1;
% flag = group==0;
% flag = true(1,length(group));

corr_type = 'pearson'; %'spearman';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


data1 = cellfun(@(x) x(flag), neural_metric{con_i}, 'uni', 0);

% data2 = cellfun(@(x) x(flag), behav_metric_list, 'uni', 0);
data2 = cellfun(@(x) x, behav_metric_list, 'uni', 0);
% data2 = {cov};

figure('position', [50 50 560 910]);
jh_corr_table(data1, data2, ...
            roi_name_list, cellfun(@(x) replace(x,'_',' '),behav_metric_name_list, 'uni', 0),  0.05, corr_type);
figure('position', [600 50 560 910]);
jh_corr_table(data1, data2, ...
            roi_name_list, cellfun(@(x) replace(x,'_',' '),behav_metric_name_list, 'uni', 0),  0.1, corr_type);



% figure;
% jh_regress(data2{3}, data1{0+5})


%% single bar

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
con_i = 1;

% flag = group==0;
flag = true(1,length(group));
flag = group == 1;

roi_list = 1:3; fig_size = [12,8]; % big ROIs
roi_list = 5:7; fig_size = [12,8]; % subregions
roi_list = 8:11; fig_size = [14,8]; % HBT

roi_list = [roi_list, roi_list+11, roi_list+22];

stat_type = 'nonparam';

color = COLOR_SPATIAL;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data = cellfun(@(x) x(flag), neural_metric{con_i}, 'uni', 0);
data = data(roi_list);

figure;
fig = jh_bar(data,'color',color, 'DrawPoint',false, 'drawstats',true, 'stattype', stat_type);
xticklabels(roi_name_list(roi_list))
xticks(find(~cellfun(@isempty, data)));
xlim([0,length(data)+1])
xline([length(data)/3+0.5, length(data)/3*2+0.5], '--r', 'linewidth',1.5)
% jh_set_fig('size',fig_size)


temp_roi_name_list = roi_name_list(roi_list);
fprintf('\n\n')
for data_i = 1:length(data)
    [~,p,~,stats] = ttest(data{data_i});
    fprintf('%-20s p = %.3f, t(%d) = %.3f\n', temp_roi_name_list{data_i},p, stats.df, stats.tstat)
end
fprintf('\n\n')

%% group comparison
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
con_i = 1;

flag1 = group==1; 
flag2 = group==0;

roi_list = 1:3; fig_size = [12,8]; % big ROIs
% roi_list = 5:7; fig_size = [12,8]; % subregions
% roi_list = 8:11; fig_size = [12,8]; % HBT

roi_list = [roi_list, roi_list+11, roi_list+22];
% roi_list = [roi_list, roi_list+11, roi_list+22];

stat_type = 'nonparam';

color = COLOR_SPATIAL;

color = repmat( {color, jh_color_modify_hsv(color, 'saturation',0.25)}, 1, length(roi_list));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data1 = cellfun(@(x) x(flag1), neural_metric{con_i}, 'uni', 0);
data2 = cellfun(@(x) x(flag2), neural_metric{con_i}, 'uni', 0);
data = [data1(roi_list);data2(roi_list)];
data = reshape(data, 1, []);

fig_size = [fig_size(1)*1.5, fig_size(2)];

figure;
fig = jh_bar(data,'color',color, 'DrawPoint',false, 'drawstats',true, 'stattype',stat_type);
% fig = jh_bar(data,'color',color);
xticklabels(roi_name_list(roi_list))
xticks(1:2:length(data));
xlim([0,length(data)+1])
xline([length(data)/3+0.5, length(data)/3*2+0.5], '--r', 'linewidth',1.5)
jh_set_fig('size',fig_size)


temp_roi_name_list = roi_name_list(roi_list);
fprintf('\n\n')
for data_i = 1:length(data)
%     [~,p,~,stats] = ttest(data{data_i});
    [p,~,stats] = signrank(data{data_i}, 0, 'method','approximate');
    fprintf('%-20s p = %.3f, z = %.3f, n = %d', temp_roi_name_list{ceil(data_i/2)}, ...
                                                  p, stats.zval, sum(~isnan(data{data_i})));
    if p < 0.05; fprintf('\t*\n'); else; fprintf('\n'); end
end
fprintf('\n\n')

temp_roi_name_list = roi_name_list(roi_list);
fprintf('\n\n')
for data_i = 1:(length(data)/2)
    temp1 = data{(data_i-1)*2+1};
    temp2 = data{(data_i-1)*2+2};
    [p,~,stats] = ranksum(temp1, temp2, 'method','approximate');
    fprintf('%-20s p = %.3f, z = %.3f', temp_roi_name_list{data_i}, ...
                                                  p, stats.zval );
    if p < 0.05; fprintf('\t*\n'); else; fprintf('\n'); end
end
fprintf('\n\n')

%% group comparison 2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
con_i = 7;

flag1 = group==1; 
flag2 = group==0;

roi_list = 1:3; fig_size = [12,8]; % big ROIs
roi_list = 5:7; fig_size = [12,8]; % subregions
roi_list = 8:11; fig_size = [14,8]; % HBT

roi_list = [roi_list, roi_list+11, roi_list+22];

stat_type = 'nonparam';

color = COLOR_SPATIAL;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data = cellfun(@(x) nanmean(x(flag1))-nanmean(x(flag2)), neural_metric{con_i}, 'uni', 0);
data = data(roi_list);

if strcmp(stat_type,'param')
    [~,p] = cellfun(@(x) ttest2(x(flag1), x(flag2)), neural_metric{con_i} );
else
    p = cellfun(@(x) ranksum(x(flag1), x(flag2)), neural_metric{con_i} );
end
p = p(roi_list);

figure;
fig = jh_bar(data,'color',color, 'DrawPoint',false, 'drawstats',true,'stats',p);
xticklabels(roi_name_list(roi_list))
xticks(find(~cellfun(@isempty, data)));
xlim([0,length(data)+1])
xline([length(data)/3+0.5, length(data)/3*2+0.5], '--r', 'linewidth',1.5)
jh_set_fig('size',fig_size)


%% correlation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
con_i = 1;

% data_behav = spatial_metric_training; data_behav_cov = spatial_metric_sess_wise{1}; is_residual = true;
% data_behav = spatial_metric_sess_wise{8}; is_residual = false;

data_behav = behav_metric_list{4}; is_residual = false;


roi_list = 1:3; fig_size = [12,8]; % big ROIs
roi_list = 5:7; fig_size = [12,8]; % subregions
% roi_list = 8:11; fig_size = [12,8]; % HBT

opt_side = 3;   % 1 2 3  LRB

roi_list = roi_list + (opt_side-1)*11;

flag = group==1;


color = COLOR_SPATIAL;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if is_residual
%     data_behav = func_get_residual(data_behav_cov(flag), data_behav(flag));
% else
%     data_behav = data_behav(flag);
% end
if is_residual
    data_behav = func_get_residual(data_behav_cov, data_behav);
else
    data_behav = data_behav;
end


for roi_i = 1:length(roi_list)

    roi = roi_list(roi_i);
    data_neural = neural_metric{con_i}{roi}(flag);

    data1 = data_behav;
    data2 = data_neural;

    figure;
    [r1,p1] = jh_regress(data1, data2,'on', ...
                        'MarkerColor', color,'markeralpha',1, ...
                        'ShadeColor', color, 'ShadeAlpha', .15 );
    [r2,p2] = jh_regress(data1, data2,'off','type','spearman');
    jh_set_fig()
    
    title_str1 = sprintf('Pearson r = %.3f (p = %.3f)',r1,p1);
    if p1 < 0.05; title_str1 = [title_str1,'*']; end
    title_str2 = sprintf('Spearman r = %.3f (p = %.3f)',r2,p2);
    if p2 < 0.05; title_str2 = [title_str2,'*']; end
    title_str = {title_str1, title_str2};
    fprintf('\n%s\n',roi_name_list{roi})
    fprintf('%s\n', title_str{:})

end


%% correlation - multiple

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
con_i = 4;

% data_behav = spatial_metric_training; data_behav_cov = spatial_metric_sess_wise{1}; is_residual = true;
% data_behav = spatial_metric_sess_wise{8}; is_residual = false;
data_behav = behav_metric_list{4}; is_residual = false;
% data_behav = behav_metric_list{4}; is_residual = true;


roi_list = 1:3; fig_size = [16,18]; % big ROIs
roi_list = 5:7; fig_size = [16,18]; % subregions
% roi_list = 8:11; fig_size = [20,18]; % HBT

% roi_list = [1,2,6]; fig_size = [16,18]; % ctx

roi_list = [roi_list, roi_list+11, roi_list+22];

flag = group==1;
% flag = group==0;


color = COLOR_SPATIAL;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if is_residual
    data_behav = func_get_residual(cov, data_behav);
    % data_behav = func_get_residual(data_behav_cov(flag), data_behav(flag));
else
    data_behav = data_behav;
end

    figure;
for roi_i = 1:length(roi_list)

    subplot(3,length(roi_list)/3,roi_i)

    roi = roi_list(roi_i);
    data_neural = neural_metric{con_i}{roi}(flag);

    data1 = data_behav;
    data2 = data_neural;

%     [r1,p1] = jh_regress(data1, data2,'on', ...
%                         'MarkerColor', color,'markeralpha',1, ...
%                         'ShadeColor', color, 'ShadeAlpha', .15 );
    [r2,p2] = jh_regress(data1, data2,'off','type','spearman');
    [r2,p2] = jh_regress(data1, data2,'off','type','spearman');

    if p2 >= 0.1
        [r1,p1] = jh_regress(data1, data2,'on', ...
                            'MarkerColor', color,'markeralpha',1, ...
                            'ShadeColor', color, 'ShadeAlpha', .15 );
    elseif p2 >=0.05
        [r1,p1] = jh_regress(data1, data2,'on', ...
                            'MarkerColor', color,'markeralpha',1, ...
                            'ShadeColor', color, 'ShadeAlpha', .15, 'linecolor', [0 1 0] );
    elseif p2 >=0.01
        [r1,p1] = jh_regress(data1, data2,'on', ...
                            'MarkerColor', color,'markeralpha',1, ...
                            'ShadeColor', color, 'ShadeAlpha', .15, 'linecolor',[0 1 0] );
    else
        [r1,p1] = jh_regress(data1, data2,'on', ...
                            'MarkerColor', color,'markeralpha',1, ...
                            'ShadeColor', color, 'ShadeAlpha', .15, 'linecolor', [0 1 0] );
    end

    
    title_str1 = sprintf('r = %.3f (p = %.3f)',r1,p1);
    if p1 < 0.05; title_str1 = [title_str1,'*']; end
    title_str2 = sprintf('rho = %.3f (p = %.3f)',r2,p2);
    if p2 < 0.05; title_str2 = [title_str2,'*']; end
    title_str = {roi_name_list{roi},title_str1, title_str2};
    title_str = {roi_name_list{roi}, title_str2};
    title(title_str, 'interpreter','tex','fontweight','bold')

end

    jh_set_fig('size', fig_size)






%% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data1 = spatial_metric_sess_wise{1}; data2 = spatial_metric_training; color = COLOR_SPATIAL;
data1 = spatial_metric_sess_wise{8}; data2 = spatial_metric_training; color = COLOR_SPATIAL;

% data1 = spatial_metric_sess_wise{8}; data2 = training_cov_behav_em_list{1}; color = COLOR_SPATIAL;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
[r1,p1] = jh_regress(data1, data2,'on', ...
                    'MarkerColor', color,'markeralpha',1, ...
                    'ShadeColor', color, 'ShadeAlpha', .15 );
[r2,p2] = jh_regress(data1, data2,'off','type','spearman');
% lim = [min(data1), max(data1)];
% lim = [lim(1)-diff(lim)*0.05, lim(2)+diff(lim)*0.05];
% xlim(lim);
% xlabel(replace(label1,'_',' '));
% ylabel(replace(label2,'_',' '));

jh_set_fig()

title_str1 = sprintf('Pearson r = %.3f (p = %.3f)',r1,p1);
if p1 < 0.05; title_str1 = [title_str1,'*']; end
title_str2 = sprintf('Spearman r = %.3f (p = %.3f)',r2,p2);
if p2 < 0.05; title_str2 = [title_str2,'*']; end
title_str = {title_str1, title_str2};
fprintf('%s\n', title_str{:})

%%



