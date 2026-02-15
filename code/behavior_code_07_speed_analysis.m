clear all
rng(918)

%% Set directory
data_dir = '../data';
data_name = 'sample_data.mat';
model_info = 'multi_node_model_info.mat';
chunking_info = 'chunking_info.mat';
%% Load data
load(fullfile(data_dir,data_name)) % table_all
load(fullfile(data_dir,model_info)) % node_centers, node_sizes
load(fullfile(data_dir,chunking_info)) 
addpath('visualization_toolbox')

sbj_list = unique(table_all.sbj)';
num_of_sbj = length(sbj_list);

correct_err = table_all.err;
correct_err(table_all.score == 0) = NaN;
incorrect_err = table_all.err;
incorrect_err(table_all.score>0) = NaN;


%% all averaged

ret_speed_avg = cellfun(@(x) nanmean(x),  table_all.ret_speed);
chunking_idx_use = chunking_idx;

for sbj_i = 1:num_of_sbj
    sbj_tag = (table_all.sbj == sbj_i);
    
    tag1 =sbj_tag&(proto_resp_tag');    
    tag2 =sbj_tag&(~proto_resp_tag');   

    ret_speed_sbj_on(sbj_i) = nanmean(ret_speed_avg(tag1));
    ret_speed_sbj_off(sbj_i) = nanmean(ret_speed_avg(tag2));
    
end

diffs = ret_speed_sbj_on -ret_speed_sbj_off;
[h,p,ci,stats] = ttest(diffs); d = stats.tstat / sqrt(stats.df+1);
fprintf('mean = %.3f, STD = %.3f\n', nanmean(diffs), std(diffs,'omitnan'))
fprintf('p = %.3f, t(%d) = %.3f, CI = [%.3f %.3f], Cohen''s d = %.3f\n',p, stats.df, stats.tstat, ci(1),ci(2), d)

% correlation with chunking index
figure('position',[800 800 300 300]);
[r,p] = jh_regress(chunking_idx_use,diffs);
set(gca,'LineWidth',1.2, 'FontWeight', 'normal','FontSize',15);
hold on
title(sprintf('r = %.3f, P = %.3f', r, p),'FontName','Helvetica','FontSize',17, 'FontWeight','normal')


%% Figure 3g
ret_speed_all = [];
part_len = 20;
for part_i = 1:part_len
    ret_speed_all(:,part_i) = cellfun(@(x) mean(x(round(length(x)/part_len*(part_i-1))+1:round(length(x)/part_len*part_i))), table_all.ret_speed);
    
end

tag_part = [1:part_len]
tag1_ret_speed = zeros(num_of_sbj, length(tag_part)); 
tag2_ret_speed = zeros(num_of_sbj, length(tag_part));
figure;
hold on
for part_i = tag_part
    for sbj_i = 1:num_of_sbj
        sbj_tag = table_all.sbj == sbj_i;
        
        tag1 =sbj_tag&(proto_resp_tag');    
        tag2 =sbj_tag&(~proto_resp_tag');   
    
        in_data1 = num2cell(ret_speed_all(tag1,part_i),1);
        [avg, err] = jh_mean_err(in_data1);
        tag1_ret_speed(sbj_i,part_i) = avg;

        in_data2 = num2cell(ret_speed_all(tag2,part_i),1);
        [avg, err] = jh_mean_err(in_data2);
        tag2_ret_speed(sbj_i,part_i) = avg;
    end
end

[h,p] = ttest(tag1_ret_speed - tag2_ret_speed);
in_data = num2cell(tag1_ret_speed - tag2_ret_speed,1);
[avg, err] = jh_mean_err(in_data);

figure;
hold on
shadedErrorBar(tag_part, avg, err)
yline(0, 'LineWidth',1.2)
scatter(find((p<0.1)&(p>0.05)),repmat(0.011,1,sum((p<0.1)&(p>0.05))),'SizeData',15)
scatter(find((p<0.05)&(p>0.01)),repmat(0.011,1,sum((p<0.05)&(p>0.01))),'SizeData',15)

