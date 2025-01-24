rng(918)

%% Set directory
data_dir = '../NC_data';
data_name = 'sample_data.mat';
node_model  = 'multi_node_model_info.mat';
%% Load data
load(fullfile(data_dir,data_name)) % table_all
load(fullfile(data_dir,node_model)) % node_centers, node_sizes
addpath('visualization_toolbox')

sbj_list = unique(table_all.sbj)';
num_of_sbj = length(sbj_list);

%% Calculate chunking index based on each individual's best fitting 9-node model

proto_coin_dist = []; proto_resp_dist = [];
proto_coin_tag = []; proto_resp_tag = [];
chunking_idx = [];

for sbj_i = 1:num_of_sbj
    sbj_i
    sbj_tag = table_all.sbj == sbj_i;

    % individual model info
    curr_center = node_centers{sbj_i}; 
    curr_size = node_sizes{sbj_i}; 

    [proto_coin_dist(sbj_tag), proto_coin_idx(sbj_tag)] =  min(pdist2(curr_center,[table_all.coin_x(sbj_tag), table_all.coin_y(sbj_tag)]),[],1);
    [proto_resp_dist(sbj_tag), proto_resp_idx(sbj_tag)] =  min(pdist2(curr_center,[table_all.resp_x(sbj_tag), table_all.resp_y(sbj_tag)]),[],1);

    proto_coin_tag(sbj_tag) =  proto_coin_dist(sbj_tag)<curr_size;
    proto_resp_tag(sbj_tag) =  proto_resp_dist(sbj_tag)<curr_size;

    chunking_idx(sbj_i) = (mean(proto_resp_tag(sbj_tag)) - mean(proto_coin_tag(sbj_tag)))/mean(proto_coin_tag(sbj_tag))*100;

end

save([data_dir,'/coin_resp_node_loc.mat'], 'proto_coin_tag', 'proto_resp_tag','chunking_idx')

%% Correlation between chunking index and task performance (error, accuracy, precision)
acc_thres = median(table_all.err);

err_sbj = []; 
acc_sbj = [];
prec_sbj = [];
for sbj_i = 1:num_of_sbj
    sbj_tag = table_all.sbj == sbj_i;
    err = table_all.err(sbj_tag);
    err_sbj(sbj_i) = mean(err);
    acc_sbj(sbj_i) = mean(err<acc_thres);
    prec_sbj(sbj_i) = -mean((err(err<acc_thres)-acc_thres)/acc_thres*100);
end


task_performance = prec_sbj; %err_sbj, acc_sbj

% Figure 3d
figure;
[r,p] = jh_regress(chunking_idx, task_performance)
fprintf('r = %.3f, p = %.3f\n', r, p) % or session


%% Chunking index across conditions (order, session)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% target variable for calculating chunking idx
target_variable = table_all.order; % table_all.sess;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

chunking_idx_vars = [];
data_len = length(unique(target_variable));
for sbj_i = 1:num_of_sbj
    sbj_i
    sbj_tag = table_all.sbj ==sbj_i;
   
    curr_center = node_centers{sbj_i};
    curr_size = node_sizes{sbj_i}; 

    proto_coin_dist(sbj_tag) =  min(pdist2(curr_center,[table_all.coin_x(sbj_tag), table_all.coin_y(sbj_tag)]),[],1);
    proto_resp_dist(sbj_tag) =  min(pdist2(curr_center,[table_all.resp_x(sbj_tag), table_all.resp_y(sbj_tag)]),[],1);
    proto_coin_tag(sbj_tag) =  proto_coin_dist(sbj_tag)<curr_size;
    proto_resp_tag(sbj_tag) =  proto_resp_dist(sbj_tag)<curr_size;

    for data_i = 1:data_len
        data_tag = sbj_tag & (target_variable == data_i);
        chunking_idx_vars(sbj_i,data_i) = (mean(proto_resp_tag(data_tag)) - mean(proto_coin_tag(data_tag)))/mean(proto_coin_tag(data_tag))*100;
    end
end

% Statistical testing
use_idx = chunking_idx_vars;
[h,p] = ttest(use_idx);

use_idx = num2cell(use_idx,2);
tag = 2:7;

% Fisher Z 
idx_change = cellfun(@(x) corr(tag',x(tag)','rows','pairwise'), use_idx);
[h,p,~,stats] = ttest(atanh(idx_change));
bar_data = {atanh(idx_change)};
[avg, err] = jh_mean_err(bar_data)

figure; % Figure 2c
jh_bar(avg, err, bar_data)
fprintf('Chunking index across order: p = %.3f, t(%d) = %.3f\n', p, stats.df, stats.tstat) % or session

%% Correlation between chunking index and error increase across orders (Figure 2d)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% target variable for calculating errors
target_variable = table_all.order;
use_idx = chunking_idx;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

err_data = [];
data_len = length(unique(target_variable));
for sbj_i = 1:num_of_sbj
    sbj_tag = table_all.sbj == sbj_i;
    for data_i = 1:data_len
        data_tag = (target_variable == data_i);
        err_data(sbj_i, data_i) = mean(table_all.err(sbj_tag&data_tag));
    end
end

norm_err = err_data - err_data(:,1);
chunker_list = find(use_idx>median(use_idx));

tag = 2:7;
figure;
hold on
slope = [];
for sbj_i = 1:num_of_sbj
    if ismember(sbj_i, chunker_list)
        plot_c = 'b';
    else 
        plot_c = 'k'
    end
    plot(norm_err(sbj_i,:),'LineWidth',2,'Color',plot_c)
    coefficients = polyfit(tag, err_data(sbj_i,tag), 1);
    slope(sbj_i) = corr(tag', [err_data(sbj_i,tag)]');
end
set(gca,'LineWidth',1.2, 'FontWeight', 'bold','FontSize',15);
xticks([1:8])
xticklabels({1:8})
xlim([0.5 8.5])

figure;
set(gca,'LineWidth',1.2, 'FontWeight', 'bold','FontSize',15);
[h,p] = jh_regress(use_idx,slope)
title(sprintf('p= %.3f',p),'FontName','Helvetica','FontSize',20, 'FontWeight','bold')
