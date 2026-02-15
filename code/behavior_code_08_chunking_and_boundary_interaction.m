clear all
rng(918)

%% Set directory
data_dir = '../data';
data_name = 'sample_data.mat';
chunking_info  = 'chunking_info.mat';

%% Load data
load(fullfile(data_dir,data_name)) % table_all
load(fullfile(data_dir,chunking_info)) % proto_coin_tag, proto_resp_tag, chunking_idx
addpath('visualization_toolbox')

sbj_list = unique(table_all.sbj)';
num_of_sbj = length(sbj_list);

correct_coin = table_all.err;
correct_coin(correct_coin>0.3) = NaN;
%% Boundary analysis setting
arena_size = 1; 
dist_bnd_all = [table_all.coin_x, table_all.coin_y, arena_size- table_all.coin_x, arena_size-table_all.coin_y];
dist_bnd = min(dist_bnd_all,[],2); dist_bnd_all(:,2); 
bnd_thres = median(dist_bnd); 

% tag_boundary coins
bnd_tag = dist_bnd<bnd_thres; 

%% boundary vs. inner

data1_tag = bnd_tag; data1_name = 'Boundary'; % 
data2_tag = ~bnd_tag; data2_name = 'Inner'; %

sbj_data1 = []; sbj_data2 = [];
for sbj_i = 1:num_of_sbj
    sbj_tag = table_all.sbj == sbj_list(sbj_i);
    sbj_data1(sbj_i) = nanmean(correct_coin(data1_tag&sbj_tag));
    sbj_data2(sbj_i) = nanmean(correct_coin(data2_tag&sbj_tag));
end

fprintf('%s: mean = %.3f, STD = %.3f\n', data1_name, nanmean(sbj_data1), std(sbj_data1,'omitnan'))
fprintf('%s: mean = %.3f, STD = %.3f\n', data2_name, nanmean(sbj_data2), std(sbj_data2,'omitnan'))

% data1 vs data2 
diffs = sbj_data1-sbj_data2;
[h,p,ci,stats] = ttest(diffs); d = stats.tstat / sqrt(stats.df+1);
fprintf('p = %.3f, t(%d) = %.3f, CI = [%.3f %.3f], Cohen''s d = %.3f\n',p, stats.df, stats.tstat, ci(1),ci(2), d)

% figure;
in_data = {sbj_data1, sbj_data2};
[avg err] = jh_mean_err(in_data)
figure;
jh_bar(avg ,err, in_data)
xticks([1 2])
xticklabels({data1_name, data2_name})
set(gca,'LineWidth',1.2, 'FontWeight', 'normal','FontSize',15);


% "Data" x "Order" repeated-measures ANOVA
data1 = []; data2 = [];
for sbj_i = 1:num_of_sbj
    sbj_tag = table_all.sbj == sbj_i;

    for order_i = 1:8
        order_tag = table_all.order == order_i;
        data1(sbj_i,order_i) = nanmean(correct_coin(data1_tag&sbj_tag&order_tag));
        data2(sbj_i,order_i) = nanmean(correct_coin(data2_tag&sbj_tag&order_tag));
    end
end

t = array2table([data1, data2]);
t.Properties.VariableNames = {'data1_1', 'data1_2', 'data1_3', 'data1_4', 'data1_5', 'data1_6', 'data1_7', 'data1_8','data2_1', 'data2_2', 'data2_3', 'data2_4', 'data2_5', 'data2_6', 'data2_7', 'data2_8'};
withinfactors = table(categorical({'data1','data1','data1','data1','data1','data1','data1','data1','data2','data2','data2','data2','data2','data2','data2','data2'})',categorical({'1','2','3','4','5','6','7','8','1','2','3','4','5','6','7','8'})','VariableNames',{'data','order'});

rm = fitrm(t, 'data1_1,data1_2,data1_3,data1_4,data1_5,data1_6,data1_7,data1_8,data2_1,data2_2,data2_3,data2_4,data2_5,data2_6,data2_7,data2_8 ~ 1','WithinDesign',withinfactors);
ranovatable = ranova(rm, 'WithinModel', 'data*order')

% figure;
in_data = num2cell(data2 - data1,1);
[avg err] = jh_mean_err(in_data)
figure;
jh_bar(avg ,err, in_data,'DrawStats',1)
xticks([1:8])
xticklabels({1:8})
ylabel(sprintf('%s > %s error diff.', data2_name, data1_name))
set(gca,'LineWidth',1.2, 'FontWeight', 'normal','FontSize',15);
%% Chunking node x Boundary proximity

tag1 = proto_coin_tag' & bnd_tag; mean(tag1)
tag2 = ~proto_coin_tag' & bnd_tag; mean(tag2)
tag3 = proto_coin_tag' & ~bnd_tag; mean(tag3)
tag4 = ~proto_coin_tag' & ~bnd_tag; mean(tag4)

% Figure 3a
figure; 
hold on
scatter(table_all.coin_x(tag1), table_all.coin_y(tag1),15,'filled'); 
scatter(table_all.coin_x(tag2), table_all.coin_y(tag2),15,'filled'); 
scatter(table_all.coin_x(tag3), table_all.coin_y(tag3),15,'filled'); 
scatter(table_all.coin_x(tag4), table_all.coin_y(tag4),15,'filled'); 
box on
axis equal
hold on; xlim([0 arena_size]); ylim([0 arena_size]); %ylim([arena_size/4 arena_size/4*3]); box on

err1 = []; err2 = []; err3 = []; err4= [];
for sbj_i = 1:num_of_sbj
    sbj_tag = (table_all.sbj==sbj_list(sbj_i));% & acc_tag

    err1(sbj_i) = mean(correct_coin(sbj_tag&tag1),'all','omitnan');
    err2(sbj_i) = mean(correct_coin(sbj_tag&tag2),'all','omitnan');
    err3(sbj_i) = mean(correct_coin(sbj_tag&tag3),'all','omitnan');
    err4(sbj_i) = mean(correct_coin(sbj_tag&tag4),'all','omitnan');
end


% Figure 3b
figure;
hold on;
in_data = num2cell([err1;err2;err3;err4],2); 
[avg, err] = jh_mean_err(in_data);
jh_bar(avg, err)
xticks([1:4])
xticklabels({"on bnd", "off bnd","on cen", "off cen"})
set(gca,'LineWidth',1.2, 'FontWeight', 'normal','FontSize',15);
meas_list = {'meas1', 'meas2', 'meas3', 'meas4', 'meas5', 'meas6', 'meas7', 'meas8', 'meas9', 'meas10'};
t = array2table([in_data{1}',in_data{2}',in_data{3}',in_data{4}']);
meas_num = size(t,2);
ylim([0 0.2])
t.Properties.VariableNames = meas_list(1:meas_num);
withinfactors = table(categorical({'bnd','bnd','center','center'})',categorical({'on','off','on','off'})','VariableNames',{'bnd','node'});
rm = fitrm(t,['meas1-meas',num2str(meas_num),'~1'],'WithinDesign',withinfactors);
tbl = ranova(rm, 'WithinModel', 'bnd*node') % AVNOA results 

ss_bnd = tbl.SumSq('(Intercept):bnd');
ss_error_bnd = tbl.SumSq('Error(bnd)');
eta_p_bnd = ss_bnd / (ss_bnd + ss_error_bnd);

ss_node = tbl.SumSq('(Intercept):node');
ss_error_node = tbl.SumSq('Error(node)');
eta_p_node = ss_node / (ss_node + ss_error_node);

ss_inter = tbl.SumSq('(Intercept):bnd:node');
ss_error_inter = tbl.SumSq('Error(bnd:node)');
eta_p_inter = ss_inter / (ss_inter + ss_error_inter);

fprintf('Partial Eta Squared (bnd): %.3f\n', eta_p_bnd);
fprintf('Partial Eta Squared (node): %.3f\n', eta_p_node);
fprintf('Partial Eta Squared (interaction): %.3f\n\n', eta_p_inter);

% on vs. off node error difference (inner/ boundaries)
boundary_diff = in_data{1} - in_data{2};
inner_diff = in_data{3} - in_data{4};
fprintf('inner: mean = %.3f, STD = %.3f\n', nanmean(inner_diff), std(inner_diff,'omitnan'))
fprintf('boundary: mean = %.3f, STD = %.3f\n\n', nanmean(boundary_diff), std(boundary_diff,'omitnan'))

diffs = inner_diff - boundary_diff;
[h,p,ci,stats] = ttest(diffs); d = stats.tstat / sqrt(stats.df+1);
fprintf('mean = %.3f, STD = %.3f\n', nanmean(diffs), std(diffs,'omitnan'))
fprintf('p = %.3f, t(%d) = %.3f, CI = [%.3f %.3f], Cohen''s d = %.3f\n',p, stats.df, stats.tstat, ci(1),ci(2), d)
