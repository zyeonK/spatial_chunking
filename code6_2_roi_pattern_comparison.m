%% load data
clear pattern_all con_name_list roi_name_list

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% reg_setting_name = 'task_period_basic_add_keypress'; %% will not use this
reg_setting_name = 'node_coin_binary';
reg_setting_name = 'node_target_coin_binary';
reg_setting_name = 'prepost_coin_binary'
% reg_setting_name = 'node_coin_binary_exclude_task';

ctx = '';
% ctx = '_ctx';
% ctx = '_hpc_ap';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load(sprintf('PATTERN_ORGANIZED%s_%s.mat', ctx, reg_setting_name))
% 'pattern_all','con_name_list','roi_name_list'
% pattern_all{con_i}{roi_i}{sbj_i} 
%%
group_tag = logical(group);

ps_all = {}; % {type}
% exp 8 - NaN
for roi_i = 1:length(roi_name_list)
    enc_on = pattern_all{1}{roi_i};
    enc_off = pattern_all{2}{roi_i};
    ret_on = pattern_all{3}{roi_i};
    ret_off = pattern_all{4}{roi_i};

    nan_idx = find(logical(cell2mat(cellfun(@(x) size(x,2)==1, enc_on, 'UniformOutput', false))));

    ps_all{roi_i} = [];
    
    for sbj_i = 1:length(sbj_list)
        if ismember(sbj_i,nan_idx)
            ps_all{roi_i}{sbj_i}{1} = nan;
            ps_all{roi_i}{sbj_i}{2} = nan;
        else

            sbj_enc_on = cell2mat(enc_on{sbj_i});
            sbj_enc_off = cell2mat(enc_off{sbj_i});
            sbj_ret_on = cell2mat(ret_on{sbj_i});
            sbj_ret_off = cell2mat(ret_off{sbj_i});
    
            if isempty(sbj_enc_on)
                ps_all{roi_i}{sbj_i}{1} = nan;
                ps_all{roi_i}{sbj_i}{2} = nan;
            else
                enc_on_off = corr(sbj_enc_on, sbj_enc_off,'rows','pairwise');
                ret_on_off = corr(sbj_ret_on, sbj_ret_off,'rows','pairwise');
                enc_ret_on = corr(sbj_enc_on, sbj_ret_on,'rows','pairwise');
                enc_ret_off = corr(sbj_enc_off, sbj_ret_off,'rows','pairwise');
    
                within_enc = diag(enc_on_off);
                within_ret = diag(ret_on_off);
                within_on = diag(enc_ret_on);
                within_off = diag(enc_ret_off);
    
                [i, j] = find(~eye(size(enc_on_off))); 
                across_enc = enc_on_off(sub2ind(size(enc_on_off), i, j));
                [i, j] = find(~eye(size(ret_on_off))); 
                across_ret = ret_on_off(sub2ind(size(ret_on_off), i, j));
                [i, j] = find(~eye(size(enc_ret_on))); 
                across_on = enc_ret_on(sub2ind(size(enc_ret_on), i, j));
                [i, j] = find(~eye(size(enc_ret_off))); 
                across_off = enc_ret_off(sub2ind(size(enc_ret_off), i, j));
    
                ps_all{roi_i}{sbj_i}{1} = [within_enc, within_ret, within_on, within_off];
                ps_all{roi_i}{sbj_i}{2} = [across_enc, across_ret, across_on, across_off];
            end
        end
        
    end
end


%%
for roi_i = 1:length(roi_name_list)
    within_data = cellfun(@(x) x{1}, ps_all{roi_i}, 'UniformOutput', false);
    across_data = cellfun(@(x) x{2}, ps_all{roi_i}, 'UniformOutput', false);

    for sbj_i = 1:length(sbj_list)
        if isnan(across_data{sbj_i})
            within_enc_all{roi_i}(sbj_i,:) = NaN;
            within_ret_all{roi_i}(sbj_i,:) = NaN;
            within_on_all{roi_i}(sbj_i,:) = NaN;
            within_off_all{roi_i}(sbj_i,:) = NaN;
            across_enc_all{roi_i}(sbj_i,:) = NaN;
            across_ret_all{roi_i}(sbj_i,:) = NaN;
            across_on_all{roi_i}(sbj_i,:) = NaN;
            across_off_all{roi_i}(sbj_i,:) = NaN;
        else
            within_enc_all{roi_i}(sbj_i,:) = within_data{sbj_i}(:,1);
            within_ret_all{roi_i}(sbj_i,:) = within_data{sbj_i}(:,2);
            within_on_all{roi_i}(sbj_i,:) = within_data{sbj_i}(:,3);
            within_off_all{roi_i}(sbj_i,:) = within_data{sbj_i}(:,4);
            across_enc_all{roi_i}(sbj_i,:) = across_data{sbj_i}(:,1);
            across_ret_all{roi_i}(sbj_i,:) = across_data{sbj_i}(:,2);
            across_on_all{roi_i}(sbj_i,:) = across_data{sbj_i}(:,3);
            across_off_all{roi_i}(sbj_i,:) = across_data{sbj_i}(:,4);
        end
    end
end
% stats
within_enc_all_avg = cellfun(@(x) nanmean(x,2), within_enc_all, 'UniformOutput', false);
within_ret_all_avg = cellfun(@(x) nanmean(x,2), within_ret_all, 'UniformOutput', false);
within_on_all_avg = cellfun(@(x) nanmean(x,2), within_on_all, 'UniformOutput', false);
within_off_all_avg = cellfun(@(x) nanmean(x,2), within_off_all, 'UniformOutput', false);
across_enc_all_avg = cellfun(@(x) nanmean(x,2), across_enc_all, 'UniformOutput', false);
across_ret_all_avg = cellfun(@(x) nanmean(x,2), across_ret_all, 'UniformOutput', false);
across_on_all_avg = cellfun(@(x) nanmean(x,2), across_on_all, 'UniformOutput', false);
across_off_all_avg = cellfun(@(x) nanmean(x,2), across_off_all, 'UniformOutput', false);

%%
sbj_tag = logical(group);
tag_left = 1:(length(roi_name_list)/3); tag_right = (length(roi_name_list)/3+1):(length(roi_name_list)/3*2);  tag_bi = (length(roi_name_list)/3*2+1):length(roi_name_list);
color_list = [repmat({"#FA8072"	},1,length(roi_list)/3),repmat({'#6495ED'},1,length(roi_list)/3),repmat({'#A9A9A9'},1,length(roi_list)/3)]

use_data = cell2mat(across_off_all_avg);
use_data(:,tag_bi) = (use_data(:,tag_left)+use_data(:,tag_right))/2;
use_data = num2cell(use_data,1);
use_data = cellfun(@(x) x(sbj_tag), use_data(roi_list), 'UniformOutput', false);
[avg, err] = jh_mean_err(use_data);
figure;
jh_bar(avg, err, use_data, 'Color', color_list, 'DrawPoint',true, 'drawstats',true,'PointSize',5)
xticks([1:length(roi_name_list(roi_list))])
xticklabels(roi_name_list(roi_list))
ylim([-0.2 0.2])
set(gca,'LineWidth',1.2, 'FontWeight', 'bold','FontSize',13);



%%

% within enc ~ trial-wise chunking
% within ret ~ trial-wise chunking
% within on ~ trial-wise chunking
% within off ~ trial-wise chunking


% within enc - within ret
% across enc - across ret
% within on - within off
% across on - across off

for roi_i = 1:length(roi_name_list)
    within_data = cellfun(@(x) x{1}, ps_all{roi_i}, 'UniformOutput', false);
    across_data = cellfun(@(x) x{2}, ps_all{roi_i}, 'UniformOutput', false);

    for sbj_i = 1:length(sbj_list)
        if isnan(across_data{sbj_i})
            within_enc_min_ret{roi_i}(sbj_i,:) = NaN;
            within_on_min_off{roi_i}(sbj_i,:) = NaN;
            across_enc_min_ret{roi_i}(sbj_i,:) = NaN;
            across_on_min_off{roi_i}(sbj_i,:) = NaN;
        else
            curr_sbj = within_data{sbj_i};
            within_enc_min_ret{roi_i}(sbj_i,:) = curr_sbj(:,1) - curr_sbj(:,2);
            within_on_min_off{roi_i}(sbj_i,:) = curr_sbj(:,3) - curr_sbj(:,4);
    
            curr_sbj = across_data{sbj_i};
            
            across_enc_min_ret{roi_i}(sbj_i,:) = curr_sbj(:,1) - curr_sbj(:,2);
            across_on_min_off{roi_i}(sbj_i,:) = curr_sbj(:,3) - curr_sbj(:,4);
        end
    end
end

% stats
within_enc_ret_avg = cellfun(@(x) nanmean(x,2), within_enc_min_ret, 'UniformOutput', false);
within_on_off_avg = cellfun(@(x) nanmean(x,2), within_on_min_off, 'UniformOutput', false);
across_enc_ret_avg = cellfun(@(x) nanmean(x,2), across_enc_min_ret, 'UniformOutput', false);
across_on_off_avg = cellfun(@(x) nanmean(x,2), across_on_min_off, 'UniformOutput', false);

%%

sbj_tag = logical(group);
tag_left = 1:(length(roi_name_list)/3); tag_right = (length(roi_name_list)/3+1):(length(roi_name_list)/3*2);  tag_bi = (length(roi_name_list)/3*2+1):length(roi_name_list);
color_list = [repmat({"#FA8072"	},1,length(roi_list)/3),repmat({'#6495ED'},1,length(roi_list)/3),repmat({'#A9A9A9'},1,length(roi_list)/3)]

use_data = cell2mat(across_on_off_avg);
use_data(:,tag_bi) = (use_data(:,tag_left)+use_data(:,tag_right))/2;
use_data = num2cell(use_data,1);
use_data = cellfun(@(x) x(sbj_tag), use_data(roi_list), 'UniformOutput', false);
[avg, err] = jh_mean_err(use_data);
figure;
jh_bar(avg, err, use_data, 'Color', color_list, 'DrawPoint',true, 'drawstats',true,'PointSize',5)
xticks([1:length(roi_name_list(roi_list))])
xticklabels(roi_name_list(roi_list))
ylim([-0.2 0.2])
set(gca,'LineWidth',1.2, 'FontWeight', 'bold','FontSize',13);
