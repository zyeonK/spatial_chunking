rng(918)

%% Set directory
data_dir = '../NC_data';
data_name = 'sample_data.mat';
coin_resp_node_loc  = 'coin_resp_node_loc.mat';

%% Load data
load(fullfile(data_dir,data_name)) % table_all
load(fullfile(data_dir,coin_resp_node_loc)) % proto_coin_tag, proto_resp_tag, chunking_idx
addpath('visualization_toolbox')

sbj_list = unique(table_all.sbj)';
num_of_sbj = length(sbj_list);

% accuracy threshold: median cut
acc_thres = median(table_all.err); 

%% Boundary analysis setting
arena_size = 1; 
dist_bnd_all = [table_all.coin_x, table_all.coin_y, arena_size- table_all.coin_x, arena_size-table_all.coin_y];
dist_bnd = min(dist_bnd_all,[],2); dist_bnd_all(:,2); 
bnd_thres = median(dist_bnd); 

% tag_boundary coins
bnd_tag = dist_bnd<bnd_thres; 

%% Boundary vs center memory difference

correct_coin = table_all.err<acc_thres; 

sbj_bnd = []; sbj_center = [];
for sbj_i = 1:num_of_sbj
    sbj_tag = table_all.sbj == sbj_i;

    sbj_bnd(sbj_i) = mean(correct_coin(bnd_tag&sbj_tag))*100;
    sbj_center(sbj_i) = mean(correct_coin(~bnd_tag&sbj_tag))*100;
end

% Boundary vs. inner
[h,p,ci, stats] = ttest(sbj_bnd, sbj_center)
fprintf('Boundary vs. innrr t-test: p = %.3f, t(%d) = %.3f\n', p, stats.df, stats.tstat)

% "Boundary" x "Order" repeated-measures ANOVA
data_bnd = []; data_center = [];
for sbj_i = 1:num_of_sbj
    sbj_tag = table_all.sbj == sbj_i;

    for order_i = 1:8
        data_tag = table_all.order == order_i;
        data_bnd(sbj_i,order_i) = mean(correct_coin(bnd_tag&sbj_tag&data_tag))*100;
        data_center(sbj_i,order_i) = mean(correct_coin(~bnd_tag&sbj_tag&data_tag))*100;
    end
end

t = array2table([data_bnd, data_center]);
t.Properties.VariableNames = {'bnd_1', 'bnd_2', 'bnd_3', 'bnd_4', 'bnd_5', 'bnd_6', 'bnd_7', 'bnd_8','cen_1', 'cen_2', 'cen_3', 'cen_4', 'cen_5', 'cen_6', 'cen_7', 'cen_8'};
withinfactors = table(categorical({'bnd','bnd','bnd','bnd','bnd','bnd','bnd','bnd','center','center','center','center','center','center','center','center'})',categorical({'1','2', '3','4','5','6','7','8','1','2','3','4','5','6','7','8'})','VariableNames',{'boundary','order'});

rm = fitrm(t, 'bnd_1,bnd_2,bnd_3,bnd_4,bnd_5,bnd_6,bnd_7,bnd_8,cen_1,cen_2,cen_3,cen_4,cen_5,cen_6, cen_7, cen_8 ~ 1','WithinDesign',withinfactors);
ranovatable = ranova(rm, 'WithinModel', 'boundary*order');


%%
tag1 = proto_coin_tag' & bnd_tag; mean(tag1)
tag2 = ~proto_coin_tag' & bnd_tag; mean(tag2)
tag3 = proto_coin_tag' & ~bnd_tag; mean(tag3)
tag4 = ~proto_coin_tag' & ~bnd_tag; mean(tag4)

acc_tag = table_all.err<acc_thres;

figure; 
hold on
scatter(table_all.coin_x(tag1), table_all.coin_y(tag1),15,'filled'); 
scatter(table_all.coin_x(tag2), table_all.coin_y(tag2),15,'filled'); 
scatter(table_all.coin_x(tag3), table_all.coin_y(tag3),15,'filled'); 
scatter(table_all.coin_x(tag4), table_all.coin_y(tag4),15,'filled'); 
box on
axis equal
hold on; xlim([0 arena_size]); ylim([0 arena_size]); %ylim([arena_size/4 arena_size/4*3]); box on

acc1 = []; acc2 = []; acc3 = []; acc4= [];
prec1 = []; prec2 = []; prec3 = []; prec4= [];
for sbj_i = 1:num_of_sbj
    sbj_tag = (table_all.sbj==sbj_i);% & acc_tag

    %acc
    acc1(sbj_i) = mean(table_all.err(sbj_tag&tag1)<acc_thres,'all','omitnan');
    acc2(sbj_i) = mean(table_all.err(sbj_tag&tag2)<acc_thres,'all','omitnan');
    acc3(sbj_i) = mean(table_all.err(sbj_tag&tag3)<acc_thres,'all','omitnan');
    acc4(sbj_i) = mean(table_all.err(sbj_tag&tag4)<acc_thres,'all','omitnan');

    %fine modulation
    prec1(sbj_i) = -(mean(table_all.err(sbj_tag&tag1&acc_tag)-acc_thres)/acc_thres*100);
    prec2(sbj_i) = -(mean(table_all.err(sbj_tag&tag2&acc_tag)-acc_thres)/acc_thres*100);
    prec3(sbj_i) = -(mean(table_all.err(sbj_tag&tag3&acc_tag)-acc_thres)/acc_thres*100);
    prec4(sbj_i) = -(mean(table_all.err(sbj_tag&tag4&acc_tag)-acc_thres)/acc_thres*100);
end

figure;
hold on
in_data = num2cell([acc1;acc2;acc3;acc4],2) % num2cell([prec1;prec2;prec3;prec4],2)
% in_data = 
[avg, err] = jh_mean_err(in_data);
jh_bar(avg, err)
xticks([1:4])
xticklabels({"on bnd", "off bnd","on cen", "off cen"})

meas_list = {'meas1', 'meas2', 'meas3', 'meas4', 'meas5', 'meas6', 'meas7', 'meas8', 'meas9', 'meas10'};
t = array2table([in_data{1}',in_data{2}',in_data{3}',in_data{4}']);
meas_num = size(t,2);

t.Properties.VariableNames = meas_list(1:meas_num);
withinfactors = table(categorical({'bnd','bnd','center','center'})',categorical({'on','off','on','off'})','VariableNames',{'bnd','node'});

rm = fitrm(t,['meas1-meas',num2str(meas_num),'~1'],'WithinDesign',withinfactors);
tbl = ranova(rm, 'WithinModel', 'bnd*node');% AVNOA results 
