% Set path
data_dir = '../data';
reg_dir = '../data/neural_regressors';
reg_task_dir = '../data/neural_task_period_info';
reg_movement_dir = '../data/neural_movement';

% load data
load(fullfile(data_dir, 'node_indiv.mat')) % 'node_coord_all','node_sizes'
node_center = node_coord_all;
node_size = cell2mat(node_sizes);

load(fullfile(data_dir, 'fMRI_sess_behav_and_marker.mat')) % 'target_coin','response','enc_time','ret_time'

% parameters
sbj_list = split(num2str(1:17));
TR = 2;
 

%% Make regressor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
reg_name = 'node_four_type';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

out_path = fullfile(reg_dir,reg_name);
out_path1 = sprintf('%s_conv',out_path);
out_path2 = sprintf('%s_original',out_path);

if ~exist(out_path1,'dir')
    mkdir(out_path1)
end
if ~exist(out_path2,'dir')
    mkdir(out_path2)
end


for sbj_i = 1:length(sbj_list)

	% load data
    sbj_name = sbj_list{sbj_i};

    sbj_target_coin = target_coin{sbj_i};
    sbj_response = response{sbj_i};
    sbj_enc_time = enc_time{sbj_i};
    sbj_ret_time = ret_time{sbj_i};

    load(sprintf('%s/%d.mat',reg_task_dir,sbj_i)); % names onsets durations
    reg_movement = readmatrix(sprintf('%s/%d.txt',reg_movement_dir,sbj_i));
    reg_template = zeros(size(reg_movement,1),1);


    % get regressor
    regressors_enc = {};
    regressors_ret = {};

    for trial_i = 1:length(sbj_target_coin)
        enc_target = sbj_target_coin{trial_i};
        ret_target = sbj_response{trial_i};
        enc_time_trial = sbj_enc_time{trial_i};
        ret_time_trial = sbj_ret_time{trial_i};
        dist_err = sqrt(sum((enc_target - ret_target).^2,2));
        
        % type each target coin/ response
        enc_bnd = []; ret_bnd = [];
        on_on = []; on_off = [];
        off_on = []; off_off = [];
        enc_on = []; ret_on = [];

        dist_bnd_all = [enc_target(:,1), enc_target(:,2), 1 - enc_target(:,1), 1 - enc_target(:,2)];
        dist_bnd = min(dist_bnd_all,[],2); 
        enc_bnd = dist_bnd;

        dist_bnd_all = [ret_target(:,1), ret_target(:,2), 1 - ret_target(:,1), 1 - ret_target(:,2)];
        dist_bnd = min(dist_bnd_all,[],2); 
        ret_bnd = dist_bnd;
        
        for coin_i = 1:length(ret_target)
            temp_dist = cellfun(@(x) norm(enc_target(coin_i,:) - x), node_center{sbj_i});
            temp_dist = min(temp_dist);
            enc_on(coin_i) = temp_dist < node_size(sbj_i);

            temp_dist = cellfun(@(x) norm(ret_target(coin_i,:) - x), node_center{sbj_i});
            temp_dist = min(temp_dist);
            ret_on(coin_i) = temp_dist < node_size(sbj_i);
        end

        on_on = enc_on&ret_on;
        on_off = enc_on&(~ret_on);
        off_on = (~enc_on)&ret_on;
        off_off = (~enc_on)&(~ret_on);

        on_on_num{sbj_i}(trial_i) = sum(on_on);
        on_off_num{sbj_i}(trial_i) = sum(on_off);
        off_on_num{sbj_i}(trial_i) = sum(off_on);
        off_off_num{sbj_i}(trial_i) = sum(off_off);
        

        enc_onset = onsets{2+8*(trial_i-1)};
        ret_onset = onsets{6+8*(trial_i-1)};
        
        % enc 1 (on-on)
        temp = reg_template;
        for coin_i = find(on_on(:))'
            idx = [enc_time_trial{coin_i}(1), enc_time_trial{coin_i}(end)] + enc_onset;
            idx = round(idx/TR);
            if idx(1)<1; idx(1) = 1; end
            if idx(end)>length(temp); idx(end) = length(temp); end
            temp(idx(1):idx(end)) = 1;
        end
        regressors_enc{end+1} = temp;
            
        
        % enc 0 (off-on)
        temp = reg_template;
        for coin_i = find(off_on(:))'
            idx = [enc_time_trial{coin_i}(1), enc_time_trial{coin_i}(end)] + enc_onset;
            idx = round(idx/TR);
            if idx(1)<1; idx(1) = 1; end
            if idx(end)>length(temp); idx(end) = length(temp); end
            temp(idx(1):idx(end)) = 1;
        end
        regressors_enc{end+1} = temp;

        % enc 3 (on-off)
        temp = reg_template;
        for coin_i = find(on_off(:))'
            idx = [enc_time_trial{coin_i}(1), enc_time_trial{coin_i}(end)] + enc_onset;
            idx = round(idx/TR);
            if idx(1)<1; idx(1) = 1; end
            if idx(end)>length(temp); idx(end) = length(temp); end
            temp(idx(1):idx(end)) = 1;
        end
        regressors_enc{end+1} = temp;
            
        
        % enc 2 (off_dff)
        temp = reg_template;
        for coin_i = find(off_off(:))'
            idx = [enc_time_trial{coin_i}(1), enc_time_trial{coin_i}(end)] + enc_onset;
            idx = round(idx/TR);
            if idx(1)<1; idx(1) = 1; end
            if idx(end)>length(temp); idx(end) = length(temp); end
            temp(idx(1):idx(end)) = 1;
        end
        regressors_enc{end+1} = temp;

        % enc 2 (boundary)
        temp = reg_template;
        for coin_i = 1:8
            idx = [enc_time_trial{coin_i}(1), enc_time_trial{coin_i}(end)] + enc_onset;
            idx = round(idx/TR);
            if idx(1)<1; idx(1) = 1; end
            if idx(end)>length(temp); idx(end) = length(temp); end
            temp(idx(1):idx(end)) = enc_bnd(coin_i);
        end
        regressors_enc{end+1} = temp;

        % enc 2 (distance error)
        temp = reg_template;
        for coin_i = 1:8
            idx = [enc_time_trial{coin_i}(1), enc_time_trial{coin_i}(end)] + enc_onset;
            idx = round(idx/TR);
            if idx(1)<1; idx(1) = 1; end
            if idx(end)>length(temp); idx(end) = length(temp); end
            temp(idx(1):idx(end)) = dist_err(coin_i);
        end
        regressors_enc{end+1} = temp;

        % ret 1 (on-on)
        temp = reg_template;
        for coin_i = find(on_on(:))'
            idx = [ret_time_trial{coin_i}(1), ret_time_trial{coin_i}(end)] + ret_onset;
            idx = round(idx/TR);
            if idx(1)<1; idx(1) = 1; end
            if idx(end)>length(temp); idx(end) = length(temp); end
            temp(idx(1):idx(end)) = 1;
        end
        regressors_ret{end+1} = temp;

        % ret 0 (off-on)
        temp = reg_template;
        for coin_i = find(off_on(:))'
            idx = [ret_time_trial{coin_i}(1), ret_time_trial{coin_i}(end)] + ret_onset;
            idx = round(idx/TR);
            if idx(1)<1; idx(1) = 1; end
            if idx(end)>length(temp); idx(end) = length(temp); end
            temp(idx(1):idx(end)) = 1;
        end
        regressors_ret{end+1} = temp;

        % ret 3 (on-off)
        temp = reg_template;
        for coin_i = find(on_off(:))'
            idx = [ret_time_trial{coin_i}(1), ret_time_trial{coin_i}(end)] + ret_onset;
            idx = round(idx/TR);
            if idx(1)<1; idx(1) = 1; end
            if idx(end)>length(temp); idx(end) = length(temp); end
            temp(idx(1):idx(end)) = 1;
        end
        regressors_ret{end+1} = temp;

        % ret 2 (off-off)
        temp = reg_template;
        for coin_i = find(off_off(:))'
            idx = [ret_time_trial{coin_i}(1), ret_time_trial{coin_i}(end)] + ret_onset;
            idx = round(idx/TR);
            if idx(1)<1; idx(1) = 1; end
            if idx(end)>length(temp); idx(end) = length(temp); end
            temp(idx(1):idx(end)) = 1;
        end
        regressors_ret{end+1} = temp;

        % ret  (boundary)
        temp = reg_template;
        for coin_i = 1:8
            idx = [ret_time_trial{coin_i}(1), ret_time_trial{coin_i}(end)] + ret_onset;
            idx = round(idx/TR);
            if idx(1)<1; idx(1) = 1; end
            if idx(end)>length(temp); idx(end) = length(temp); end
            temp(idx(1):idx(end)) = ret_bnd(coin_i);
        end
        regressors_ret{end+1} = temp;

        % ret  (distance error)
        temp = reg_template;
        for coin_i = 1:8
            idx = [ret_time_trial{coin_i}(1), ret_time_trial{coin_i}(end)] + ret_onset;
            idx = round(idx/TR);
            if idx(1)<1; idx(1) = 1; end
            if idx(end)>length(temp); idx(end) = length(temp); end
            temp(idx(1):idx(end)) = dist_err(coin_i);
        end
        regressors_ret{end+1} = temp;

    end
    regressors_enc = cell2mat(regressors_enc);
    regressors_ret = cell2mat(regressors_ret);

    regressors_enc_final = []; regressors_ret_final = [];


    enc_reg_len = 6;
    regressors_enc_final(:,1) = sum(regressors_enc(:,[1,1+enc_reg_len,1+2*enc_reg_len]),2);
    regressors_enc_final(:,2) = sum(regressors_enc(:,[2,2+enc_reg_len,2+2*enc_reg_len]),2);
    regressors_enc_final(:,3) = sum(regressors_enc(:,[3,3+enc_reg_len,3+2*enc_reg_len]),2);
    regressors_enc_final(:,4) = sum(regressors_enc(:,[4,4+enc_reg_len,4+2*enc_reg_len]),2);
    regressors_enc_final(:,5) = sum(regressors_enc(:,[5,5+enc_reg_len,5+2*enc_reg_len]),2);
    regressors_enc_final(:,6) = sum(regressors_enc(:,[6,6+enc_reg_len,6+2*enc_reg_len]),2);

    ret_reg_len = 6;
    regressors_ret_final(:,1) = sum(regressors_ret(:,[1,1+ret_reg_len,1+2*ret_reg_len]),2);
    regressors_ret_final(:,2) = sum(regressors_ret(:,[2,2+ret_reg_len,2+2*ret_reg_len]),2);
    regressors_ret_final(:,3) = sum(regressors_ret(:,[3,3+ret_reg_len,3+2*ret_reg_len]),2);
    regressors_ret_final(:,4) = sum(regressors_ret(:,[4,4+ret_reg_len,4+2*ret_reg_len]),2);
    regressors_ret_final(:,5) = sum(regressors_ret(:,[5,5+ret_reg_len,5+2*ret_reg_len]),2);
    regressors_ret_final(:,6) = sum(regressors_ret(:,[6,6+ret_reg_len,6+2*ret_reg_len]),2);


    % convolution
    regressors_enc_conv = [];
    for i = 1:size(regressors_enc_final,2)
        temp = conv(regressors_enc_final(:,i),spm_hrf(TR));
        temp = temp(1:length(reg_template));
        regressors_enc_conv(:,i) = temp;
    end
    regressors_ret_conv = [];
    for i = 1:size(regressors_ret_final,2)
        temp = conv(regressors_ret_final(:,i),spm_hrf(TR));
        temp = temp(1:length(reg_template));
        regressors_ret_conv(:,i) = temp;
    end
    
    
    % save
    writematrix([regressors_enc_conv,regressors_ret_conv], fullfile(out_path1, sprintf('%d.txt',sbj_i)), 'delimiter','tab');
    writematrix([regressors_enc_final,regressors_ret_final], fullfile(out_path2, sprintf('%d.txt',sbj_i)), 'delimiter','tab');
    
    fprintf('\n%d\n',sbj_i);
end

% Check
on_on_tag = cellfun(@(x) sum(x)==0, on_on_num);
on_off_tag = cellfun(@(x) sum(x)==0, on_off_num);
off_on_tag = cellfun(@(x) sum(x)==0, off_on_num);
off_off_tag = cellfun(@(x) sum(x)==0, off_off_num);







