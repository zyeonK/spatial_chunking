clear; clc;
%% get source data

data_7t = '../source_data/Source_Data_7T_Fig4_5.xlsx';
data_3t = '../source_data/Source_Data_3T_Fig4_5.xlsx';

onoff_3t = table2array(readtable(data_3t,'Sheet','Fig4d_on_off_retrieval'));
onoff_3t = num2cell(onoff_3t(:,2:end),1);

onoff_7t = table2array(readtable(data_7t,'Sheet','Fig4e_on_off_retrieval'));
onoff_7t = num2cell(onoff_7t(:,2:end),1);

% node attraction
corr1_3t = table2array(readtable(data_3t,'Sheet','Fig4f_corr'));
corr1X_3t = corr1_3t(:,2);
corr1Y_3t = corr1_3t(:,3:end);

corr1_7t = table2array(readtable(data_7t,'Sheet','Fig4g_corr'));
corr1X_7t = corr1_7t(:,2);
corr1Y_7t = corr1_7t(:,3:end);

% avg error
corr2_3t = table2array(readtable(data_3t,'Sheet','Fig5a_corr'));
corr2X_3t = corr2_3t(:,2);
corr2Y_3t = corr2_3t(:,3:end);

corr2_7t = table2array(readtable(data_7t,'Sheet','Fig5b_corr'));
corr2X_7t = corr2_7t(:,2);
corr2Y_7t = corr2_7t(:,3:end);


roi_name_list_3t = {'Ant HP', 'Post HP'};
roi_name_list_7t = {'Ant DG','Post DG','Ant CA23','Post CA23','Ant CA1','Post CA1','Ant Sub','Post Sub'};

%% On/off-node response activation difference (Figure 4)

figure('Color','w','Position',[60 80 1100 760]);

subplot(1,2,1);
[avg, err] = jh_mean_err(onoff_3t);
jh_bar(avg, err, onoff_3t, 'DrawPoint',1, 'DrawStats',1)
ylabel('On > Off response (a.u.)'); title('Fig 4d  (3T)');
xticks(1:length(roi_name_list_3t)); xticklabels(roi_name_list_3t);

subplot(1,2,2);
[avg, err] = jh_mean_err(onoff_7t);
jh_bar(avg, err, onoff_7t, 'DrawPoint',1, 'DrawStats',1)
ylabel('On > Off response (a.u.)'); title('Fig 4e  (7T)');
xticks(1:length(roi_name_list_7t)); xticklabels(roi_name_list_7t);


%% Node attraction idx ~ activation (Figure 4)

figure;
sgtitle('Fig 4f (3T)','FontWeight', 'bold','FontSize',15)
for roi_i = 1:length(roi_name_list_3t)

    subplot(length(roi_name_list_3t),1,roi_i)
    
    hold on;
    [r,p] = jh_regress(corr1X_3t, corr1Y_3t(:,roi_i));
    axis square
    xlabel('Activation')
    ylabel('Node attraction idx')
    set(gca,'LineWidth',1.5, 'FontWeight', 'normal','FontSize',11);
    title(sprintf('%s\n(r = %.3f, p = %.3f)', roi_name_list_3t{roi_i},r,p),'FontSize',13)

end

figure;
roi_list = [1:2:8, 2:2:8];

sgtitle('Fig 4g (7T)','FontWeight', 'bold','FontSize',15)
for roi_i = 1:length(roi_name_list_7t)

    subplot(2,length(roi_name_list_7t)/2,roi_i)
    
    hold on;
    [r,p] = jh_regress(corr1X_7t, corr1Y_7t(:,roi_list(roi_i)));
    axis square
    xlabel('Activation')
    ylabel('Node attraction idx')
    set(gca,'LineWidth',1.5, 'FontWeight', 'normal','FontSize',11);
    title(sprintf('%s\n(r = %.3f, p = %.3f)', roi_name_list_7t{roi_i},r,p),'FontSize',13)

end


%% Node attraction idx ~ activation (Figure 5)

figure;
sgtitle('Fig 5a (3T)','FontWeight', 'bold','FontSize',15)
for roi_i = 1:length(roi_name_list_3t)

    subplot(length(roi_name_list_3t),1,roi_i)
    
    hold on;
    [r,p] = jh_regress(corr2X_3t, corr2Y_3t(:,roi_i));
    axis square
    xlabel('Activation')
    ylabel('Avg. dist. err.')
    set(gca,'LineWidth',1.5, 'FontWeight', 'normal','FontSize',11);
    title(sprintf('%s\n(r = %.3f, p = %.3f)', roi_name_list_3t{roi_i},r,p),'FontSize',13)

end

figure;
roi_list = [1:2:8, 2:2:8];

sgtitle('Fig 5b (7T)','FontWeight', 'bold','FontSize',15)
for roi_i = 1:length(roi_name_list_7t)

    subplot(2,length(roi_name_list_7t)/2,roi_i)
    
    hold on;
    [r,p] = jh_regress(corr2X_7t, corr2Y_7t(:,roi_list(roi_i)));
    axis square
    xlabel('Activation')
    ylabel('Avg. dist. err.')
    set(gca,'LineWidth',1.5, 'FontWeight', 'normal','FontSize',11);
    title(sprintf('%s\n(r = %.3f, p = %.3f)', roi_name_list_7t{roi_list(roi_i)},r,p),'FontSize',13)

end