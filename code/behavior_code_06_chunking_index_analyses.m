clear all
rng(918)

%% Set directory
data_dir = '../data';
data_name = 'sample_data.mat';
model_info = 'multi_node_model_info.mat';

%% Load data
load(fullfile(data_dir,data_name)) % table_all
load(fullfile(data_dir,model_info)) % node_centers, node_sizes
addpath('visualization_toolbox')

sbj_list = unique(table_all.sbj)';
num_of_sbj = length(sbj_list);

correct_err = table_all.err;
correct_err(table_all.score == 0) = NaN;
incorrect_err = table_all.err;
incorrect_err(table_all.score>0) = NaN;

%% Calculate chunking index based on each individual's best fitting 9-node model

proto_coin_dist = []; proto_resp_dist = [];
proto_coin_tag = []; proto_resp_tag = [];
sbj_avg_perf = []; chunking_idx = [];

for sbj_i = 1:num_of_sbj
    sbj_idx = sbj_list(sbj_i);
    sbj_tag = table_all.sbj == sbj_idx;
    sbj_avg_perf(sbj_i) = nanmean(correct_err(sbj_tag));
    curr_center = node_centers{sbj_i};
    curr_size = node_sizes{sbj_i}; 

    [proto_coin_dist(sbj_tag), proto_coin_idx(sbj_tag)] =  min(pdist2(curr_center,[table_all.coin_x(sbj_tag), table_all.coin_y(sbj_tag)]),[],1);
    [proto_resp_dist(sbj_tag), proto_resp_idx(sbj_tag)] =  min(pdist2(curr_center,[table_all.resp_x(sbj_tag), table_all.resp_y(sbj_tag)]),[],1);

    proto_coin_tag(sbj_tag) =  proto_coin_dist(sbj_tag)<curr_size;
    proto_resp_tag(sbj_tag) =  proto_resp_dist(sbj_tag)<curr_size;

    chunking_idx(sbj_i) = (mean(proto_resp_tag(sbj_tag)) - mean(proto_coin_tag(sbj_tag)))/mean(proto_coin_tag(sbj_tag))*100;
    % chunking_idx(sbj_i) = (mean(proto_coin_dist(sbj_tag))/mean(proto_resp_dist(sbj_tag))) - 1;

end

save(fullfile(data_dir,'chunking_info.mat'), 'chunking_idx','proto_coin_tag','proto_resp_tag')
%% Correlation between chunking index and task performance (error)
perf_metric = sbj_avg_perf; 
perf_name = 'avgerage error';
figure;
[r,p] = jh_regress(chunking_idx, perf_metric);
set(gca,'LineWidth',1.5, 'FontWeight', 'normal','FontSize',12);
xlabel('chunking all')
ylabel(perf_name)
title(sprintf('r = %.3f\n P = %.3f',r,p),'FontName','Helvetica','FontSize',14, 'FontWeight','normal')

%% Chunking index across conditions (order, session)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% target variable for calculating chunking idx
target_variable = table_all.order; %table_all.sess; % 
tag = 2:7; %1:7 %
xlabelname = 'order'; %'sess';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

chunking_idx_vars = []; err_data = []; 
data_len = length(unique(target_variable));
for sbj_i = 1:num_of_sbj
    sbj_idx = sbj_list(sbj_i);
    sbj_tag = table_all.sbj == sbj_idx;

    for data_i = 1:data_len
        data_tag = sbj_tag & (target_variable == data_i);

        err_data(sbj_i, data_i) = nanmean(correct_err(sbj_tag&data_tag));
        % err_data(sbj_i, data_i) = nanmean(table_all.err(sbj_tag&data_tag));
        if mean(proto_coin_tag(data_tag)) == 0
            chunking_idx_vars(sbj_i,data_i) = NaN;
            continue;
        end
        chunking_idx_vars(sbj_i,data_i) = (mean(proto_resp_tag(data_tag)) - mean(proto_coin_tag(data_tag)))/mean(proto_coin_tag(data_tag))*100;
        % chunking_idx_vars(sbj_i,data_i) = (mean(proto_coin_dist(data_tag))/mean(proto_resp_dist(data_tag))) - 1;

    end
end

% Statistical testing
use_idx = chunking_idx_vars;
use_data = num2cell(use_idx,1);
[avg err]= jh_mean_err(use_data);
figure;
jh_bar(avg,err,use_data,'DrawStats',1,'DrawPoint',0)
set(gca,'LineWidth',1.5, 'FontWeight', 'normal','FontSize',14);
ylabel('Chunking idx')
xlabel(xlabelname)
xticks([1:length(avg)]); xticklabels({1:length(avg)})

% % session effect in chunking
% test_data = chunking_idx_vars;
% target_variable = 'session';
% t = array2table(test_data);
% t.Properties.VariableNames = {'meas1', 'meas2', 'meas3', 'meas4', 'meas5', 'meas6', 'meas7'};
% Meas = table([1:7]','VariableNames',{target_variable});
% rm = fitrm(t,'meas1-meas7~1','WithinDesign',Meas);
% ranovatbl = ranova(rm);
% stats = multcompare(rm,target_variable);
% stats.pValue<0.05;
% SS_effect = ranovatbl.SumSq(1);
% SS_error  = ranovatbl.SumSq(2);
% eta_p_squared = SS_effect / (SS_effect + SS_error);
%% Error/ Chunking index change across variable

% chunking index across variable (bar graph) Figure 3d
use_idx = chunking_idx_vars;
figure;
hold on
yyaxis right
in_data = num2cell(err_data,1);
[avg,err] = jh_mean_err(in_data)
errorbar(avg,err,LineWidth=1.5)
ylim([0.05 0.15])
yyaxis left
in_data = num2cell(use_idx,1); %%%%%%
[avg, err] = jh_mean_err(in_data);
jh_bar(avg, err)
ylim([0 100])
xticks([1:data_len])
xticklabels({1:data_len})
set(gca,'LineWidth',1.5, 'FontWeight', 'normal','FontSize',15);

%%%%%% error
% Fisher Z 
use_err = num2cell(err_data,2);
err_change = atanh(cellfun(@(x) corr(tag',x(tag)','rows','pairwise'), use_err));

% err_change = cellfun(@(x) fitlm(tag',x(tag)'),use_err,UniformOutput=false);
% err_change = cellfun(@(x) x.Coefficients.Estimate(2),err_change);

[h,p,~,stats] = ttest(err_change);
bar_data = {err_change};
[avg, err] = jh_mean_err(bar_data)

figure;
jh_bar(avg, err, bar_data,'DrawStats',1)
ylabel('Corr. coeff (z)')
set(gca,'LineWidth',1.5, 'FontWeight', 'normal','FontSize',14);
xticklabels({})

fprintf('mean = %.3f, STD = %.3f\n', nanmean(err_change), std(err_change,'omitnan'))
[h,p,ci,stats] = ttest(err_change); %% loo_ps_r_all
d = stats.tstat / sqrt(stats.df+1);
fprintf('Error across %s: p = %.3f, t(%d) = %.3f, CI = [%.3f %.3f], Cohen''s d = %.3f\n', xlabelname, p, stats.df, stats.tstat, ci(1),ci(2), d)


%%%%%% chunking index
% Fisher Z 
use_idx = num2cell(chunking_idx_vars,2);
idx_change = atanh(cellfun(@(x) corr(tag',x(tag)','rows','pairwise'), use_idx));

% idx_change = cellfun(@(x) fitlm(tag',x(tag)'),use_idx,UniformOutput=false);
% idx_change = cellfun(@(x) x.Coefficients.Estimate(2),idx_change);

[h,p,~,stats] = ttest(idx_change);
bar_data = {idx_change};
[avg, err] = jh_mean_err(bar_data)

figure; % Figure 3e
jh_bar(avg, err, bar_data,'DrawStats',1)
ylabel('Corr. coeff (z)')
set(gca,'LineWidth',1.5, 'FontWeight', 'normal','FontSize',14);
xticklabels({})

fprintf('mean = %.3f, STD = %.3f\n', nanmean(idx_change), std(idx_change,'omitnan'))
[h,p,ci,stats] = ttest(idx_change); %% loo_ps_r_all
d = stats.tstat / sqrt(stats.df+1);
fprintf('chunking across %s: p = %.3f, t(%d) = %.3f, CI = [%.3f %.3f], Cohen''s d = %.3f\n', xlabelname, p, stats.df, stats.tstat, ci(1),ci(2), d)


%% Correlation between chunking index and error increase across orders (Figure 3f)

norm_err = err_data - err_data(:,1);
chunker_list = find(chunking_idx>median(chunking_idx));

figure;
hold on
slope = [];
for sbj_i = 1:num_of_sbj
    if ismember(sbj_i, chunker_list)
        plot_c = 'r';
    else 
        plot_c = 'k'
    end
    plot(norm_err(sbj_i,:),'LineWidth',2,'Color',plot_c)
    % polyfit_data = polyfit(tag, err_data(sbj_i,tag), 1);
    % slope(sbj_i) = polyfit_data(1);
    slope(sbj_i) = corr(tag', [err_data(sbj_i,tag)]');
end
set(gca,'LineWidth',1.2, 'FontWeight', 'normal','FontSize',12);
xticks([1:8])
xticklabels({1:8})
xlim([0.5 8.5])
xlabel('coin order')
ylabel('Distance error (normalized)')

figure;
[r,p] = jh_regress(chunking_idx,atanh(slope))
% [h,p] = jh_regress(use_idx,slope)
xlabel('Chunking indxe')
% ylabel('Error increase (r)')
ylabel('Error increase (z)')
set(gca,'LineWidth',1.2, 'FontWeight', 'normal','FontSize',12);
title(sprintf('R = %.3f\n p = %.3f',r, p),'FontName','Helvetica','FontSize',14, 'FontWeight','normal')

