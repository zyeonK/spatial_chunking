load('Z:\Spatial\data\Behavior\processed\training\training_table_all_final_enc_resp_sess8.mat');

table_all = last_sess_only

%%
load('Z:\Spatial\data\Behavior\processed\training\paper\table_all.mat') % table_all
num_of_sbj = 17;
%% Basic analysis
cov_err = []; cov_acc = []; cov_prec = [];
acc_thres = median(table_all.err);
for sbj_i = 1:num_of_sbj
    sbj_tag = table_all.sbj == sbj_i;
    cov_err(sbj_i) = mean(table_all.err(sbj_tag)); % err
    cov_acc(sbj_i) = mean(table_all.err(sbj_tag)<acc_thres)*100; %accuracy

    tag_corr = table_all.err<acc_thres;
    cov_prec(sbj_i) = mean((acc_thres - table_all.err(sbj_tag&tag_corr))/acc_thres)*100; %precision
end


%across variables
data_name= 'order'
curr_data = table_all.order;
data_list = unique(curr_data)';
err_all = [];
for sbj_i = 1:length(sbj_list)
    sbj_tag = table_all.sbj == sbj_i;
    for data_i = 1:length(data_list)
        err_data = table_all.err(sbj_tag & (curr_data == data_list(data_i)));

        err_all(sbj_i, data_i) = nanmean(err_data);
    end
end

tag_sbj =  1:17;

in_data = num2cell(err_all(tag_sbj,:),1);

figure;
hold on
set(gca,'LineWidth',1.2,'FontSize',15,'FontWeight','bold');
[avg, err] = jh_mean_err(in_data);
jh_bar(avg, err, in_data)
xlabel(data_name)
ylabel('Average error')
hold on
for sbj_i = tag_sbj
    plot(err_all(sbj_i,:))
end
xticks(1:length(data_list))
xticklabels(num2cell(data_list))


%%
t = array2table(err_all);
% t.Properties.VariableNames = {'meas1', 'meas2', 'meas3'};
% Meas = table([1 2 3 ]','VariableNames',{'Meas'});
t.Properties.VariableNames = {'meas1', 'meas2', 'meas3', 'meas4', 'meas5'};
Meas = table([1 2 3 4 5]','VariableNames',{'Meas'});
% t.Properties.VariableNames = {'meas1', 'meas2', 'meas3', 'meas4', 'meas5', 'meas6', 'meas7', 'meas8'};
% Meas = table([1 2 3 4 5 6 7 8]','VariableNames',{'Meas'});

rm = fitrm(t,'meas1-meas5~1','WithinDesign',Meas);
ranovatbl = ranova(rm)

stats = multcompare(rm,'Meas')
stats.pValue<0.05


%% boundary
% distance to boundary
arena_size = 1; 42.42;21.21; 
dist_bnd_all = [table_all.coin_x, table_all.coin_y, arena_size- table_all.coin_x, arena_size-table_all.coin_y]

figure; scatter(table_all.coin_x, table_all.coin_y,10,'filled'); axis equal
hold on; xlim([0 arena_size]); ylim([0 arena_size]); %ylim([arena_size/4 arena_size/4*3]); box on

dist_bnd = min(dist_bnd_all,[],2); dist_bnd_all(:,2); 
bnd_tag = dist_bnd<median(dist_bnd);
figure; scatter(table_all.coin_x(bnd_tag), table_all.coin_y(bnd_tag),10,'filled'); axis equal
hold on; xlim([0 arena_size]); ylim([0 arena_size]); %ylim([arena_size/4 arena_size/4*3]); box on


% Boundary vs center
bnd_thres = median(dist_bnd); % median
% bnd_thres = (1-sqrt(arena_size/2))/2 % arena size
mean(dist_bnd < bnd_thres)    

bnd_tag = dist_bnd<bnd_thres;

% performance_metric = table_all.err;
performance_metric = table_all.err<median(table_all.err);

% bar comparison

sbj_bnd = []; sbj_center = [];
for sbj_i = 1:num_of_sbj
    sbj_tag = table_all.sbj == sbj_i

    sbj_bnd(sbj_i) = mean(performance_metric(bnd_tag&sbj_tag))*100;
    sbj_center(sbj_i) = mean(performance_metric(~bnd_tag&sbj_tag))*100;
end
cov_bnd = sbj_bnd - sbj_center;
group_idx = repmat([1,2],num_of_sbj,1);
colorarray = [0.81,0.42,0.42; 0.16,0.25,0.36]
colorarray = [0.78,0.32,0; 0.07,0.44,0.67]
% colorarray = {'red','blue'};

figure('position',[800 800 200 300]);
in_data = {sbj_bnd, sbj_center}
[avg, err] = jh_mean_err(in_data)
% jh_bar(avg,err,in_data)
boxchart([sbj_bnd, sbj_center]','GroupByColor',group_idx(:))
set(gca,'LineWidth',1.2, 'FontWeight', 'bold','FontSize',15);
colororder(colorarray)
ylim([0 100])
xticklabels({})

[h,p,ci, tstats] = ttest(sbj_bnd, sbj_center)

%%
training_cov = [cov_err; cov_acc; cov_prec; cov_bnd];
save('training_cov_list.mat','cov_err','cov_acc','cov_prec','cov_bnd')

% sess8_cov = [cov_err; cov_acc; cov_prec; cov_bnd];
% save('sess8_cov_list.mat','cov_err','cov_acc','cov_prec','cov_bnd')


%%
% label_list = {'training err','training acc','training prec','training bnd','training chunking','sess8 err','sess8 acc','sess8 prec','sess8 bnd','sess8 chunking'}
% data_list = [training_cov; training_chunking; sess8_cov; sess8_chunking]';
[r,p] = corr(data_list);
p(r<0) = -p(r<0);
% corr(training_cov', sess8_cov')
figure;
jh_p_table(p, label_list,label_list,0.05)
set(gca, 'interpreter','tex','Fontsize',15)

% save('training_sess8_covariates.mat','label_list','data_list')

%%

err_pre_Post = [mean(err_all(:,1:4),2), mean(err_all(:,5:8),2)]