

%% parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sess_i = 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
path_behav = '../data/data_processed_sbj';

sbj_list = split(num2str(1:33));

TR = 2;

load('../data/data_fmri_quality_control/spatial.mat')
load('../data/data_fmri_quality_control/spatial_valid_sbj.mat') % spatial_valid_sbj

%% load data
data_all = {};
for sbj_i = 1:length(sbj_list)
    data = load(fullfile(path_behav,[sbj_list{sbj_i} '.mat']));
    data_all{sbj_i} = data.sbj;
end

age = cellfun(@(x) x.age, data_all);
sex = cellfun(@(x) x.sex, data_all);
group = cellfun(@(x) x.type,data_all);

num_sbj = length(data_all);

dir_reg_task = '../data/data_regressors_spatial/task_period';
dir_reg_movement = '../data/data_regressors_spatial/movement';

%% coin-level, binarized
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out_path = '../data/data_regressors_spatial/prototype_continuous'; type_flag = 'prototype';
out_path = '../data/data_regressors_spatial/boundary_continuous'; type_flag = 'boundary';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

out_path1 = sprintf('%s_conv',out_path);
out_path2 = sprintf('%s_original',out_path);

if ~exist(out_path1,'dir')
    mkdir(out_path1)
end
if ~exist(out_path2,'dir')
    mkdir(out_path2)
end

%
for sbj_i = 1:length(sbj_list)

	% load data
    sbj = data_all{sbj_i};
    sbj_name = sbj_list{sbj_i};
    trials = sbj.spatial.sess(end).trials;
	marker = sbj.spatial.sess(end).marker;

    load(sprintf('%s/%d.mat',dir_reg_task,sbj_i)); % names onsets durations
    reg_movement = readmatrix(sprintf('%s/%d.txt',dir_reg_movement,sbj_i));
    reg_template = zeros(size(reg_movement,1),1);


    loc_ret = cell2mat(trials(trial_i).ret.loc');
    

    % get regressor
    regressors_enc = {};
    regressors_ret = {};
    for trial_i = 1:length(trials)
        enc_onset = onsets{2+8*(trial_i-1)};
        ret_onset = onsets{6+8*(trial_i-1)};


        % enc 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        loc = cell2mat(trials(trial_i).enc.loc');
        time = cell2mat(trials(trial_i).enc.time) + enc_onset;
        if strcmp(type_flag, 'prototype')
        elseif strcmp(type_flag, 'boundary')
            [dist, idx_group] = min([loc, 1-loc],[],2);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        temp = reg_template;
        idx = round(time/TR);
        idx(idx<1) = 1;
        idx(idx>length(temp)) = length(temp);
        for i = unique(idx)
            temp(i) = mean(dist(idx==i));
        end
        regressors_enc{end+1} = temp;
        
        % enc 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        loc = cell2mat(trials(trial_i).ret.loc');
        time = cell2mat(trials(trial_i).ret.time) + ret_onset;
        if strcmp(type_flag, 'prototype')
        elseif strcmp(type_flag, 'boundary')
            [dist, idx_group] = min([loc, 1-loc],[],2);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        temp = reg_template;
        idx = round(time/TR);
        idx(idx<1) = 1;
        idx(idx>length(temp)) = length(temp);
        for i = unique(idx)
            temp(i) = mean(dist(idx==i));
        end
        regressors_ret{end+1} = temp;
    
    end
    
    regressors_enc = cell2mat(regressors_enc);
    regressors_ret = cell2mat(regressors_ret);

    % convolution
    regressors_enc_conv = [];
    for i = 1:size(regressors_enc,2)
        temp = conv(regressors_enc(:,i),spm_hrf(TR));
        temp = temp(1:length(reg_template));
        regressors_enc_conv(:,i) = temp;
    end
    regressors_ret_conv = [];
    for i = 1:size(regressors_ret,2)
        temp = conv(regressors_ret(:,i),spm_hrf(TR));
        temp = temp(1:length(reg_template));
        regressors_ret_conv(:,i) = temp;
    end
    
    
    % save
    writematrix([regressors_enc_conv,regressors_ret_conv], fullfile(out_path1, sprintf('%d.txt',sbj_i)), 'delimiter','tab');
    writematrix([regressors_enc,regressors_ret], fullfile(out_path2, sprintf('%d.txt',sbj_i)), 'delimiter','tab');
    

    fprintf('\n%d\n',sbj_i);
end



%% visualize

sbj_i = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
regressors = readmatrix(sprintf('../data/data_regressors_spatial/boundary_continuous_original/%d.txt',sbj_i));
regressors = readmatrix(sprintf('../data/data_regressors_spatial/boundary_continuous_conv/%d.txt',sbj_i));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% key press
figure('position',[50 50 1800 600]);
subplot(211)
data = regressors(:,1:3);
plot(data,'linewidth',1.2)
set(gca,'LineWidth', .8,'FontSize',11, 'FontWeight','bold');
box off
subplot(212)
data = regressors(:,4:end);
plot(data,'linewidth',1.2)
set(gca,'LineWidth', .8,'FontSize',11, 'FontWeight','bold');
box off


%%












