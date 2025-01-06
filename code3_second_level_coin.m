%% parameters

opt_img = 1; % norm_smooth / norm_rough / indiv_rough
opt_fieldmap = 1; % no-correction / correction

sess_i = 2;

opt_list = { [1 1] };


% % load behavior
% sbj_list = split(num2str(1:33));
% path_behav = '../data/data_processed_sbj/';
% data_all = {};
% for sbj_i = 1:length(sbj_list)
%     data = load(fullfile(path_behav,[sbj_list{sbj_ssdfsi} '.mat']));
%     data_all{sbj_i} = data.sbj;
% end
% group = cellfun(@(x) x.type,data_all);
% group(18) = nan; %%%%%%%%%%%%%%%


for opt_i = 1:length(opt_list)
    opt_img = opt_list{opt_i}(1);
    opt_fieldmap = opt_list{opt_i}(2);


%% settings
sbj_list = split(num2str(1:33));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
reg_setting_name = 'task_period_basic';
% reg_setting_name = 'task_period_basic_add_keypress'; 
% reg_setting_name = 'node_onoff_boundary_dist_each';
% reg_setting_name = 'prepost_coin_binary';
reg_setting_name = 'node_four_type_boundary_dist';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
dir_working = pwd();
path_behav = '../data/data_processed_sbj';

%
% con_name_list = {'enc1-0','enc0-1','enc1','enc0', ...
%                  'ret1-0','ret0-1','ret1','ret0', ...
%                  'enc_raw', 'ret_raw', 'enc_fix', 'ret_fix','enc_ctrl','ret_ctrl','enc_ctrl_fix','ret_ctrl_fix'};

% 
% con_eval_str_list = { 'arrayfun(@(x) x.enc1, beta_all{sbj_i}.beta, ''uni'', 0)', ...
%                       'arrayfun(@(x) x.enc0, beta_all{sbj_i}.beta, ''uni'', 0)', ...
%                       'arrayfun(@(x) x.enc1, beta_all{sbj_i}.beta_raw, ''uni'', 0)', ...
%                       'arrayfun(@(x) x.enc0, beta_all{sbj_i}.beta_raw, ''uni'', 0)', ...
%                       'arrayfun(@(x) x.ret1, beta_all{sbj_i}.beta, ''uni'', 0)', ...
%                       'arrayfun(@(x) x.ret0, beta_all{sbj_i}.beta, ''uni'', 0)', ...
%                       'arrayfun(@(x) x.ret1, beta_all{sbj_i}.beta_raw, ''uni'', 0)', ...
%                       'arrayfun(@(x) x.ret0, beta_all{sbj_i}.beta_raw, ''uni'', 0)', ...
%                       'arrayfun(@(x) x.enc, beta_all{sbj_i}.beta_raw, ''uni'', 0)', ...
%                       'arrayfun(@(x) x.ret, beta_all{sbj_i}.beta_raw, ''uni'', 0)', ...
%                       'arrayfun(@(x) x.enc, beta_all{sbj_i}.beta, ''uni'', 0)', ...
%                       'arrayfun(@(x) x.ret, beta_all{sbj_i}.beta, ''uni'', 0)', ...
%                       'arrayfun(@(x) x.enc_ctrl, beta_all{sbj_i}.beta, ''uni'', 0)', ...
%                       'arrayfun(@(x) x.ret_ctrl, beta_all{sbj_i}.beta, ''uni'', 0)', ...
%                       'arrayfun(@(x) x.enc_ctrl_fix, beta_all{sbj_i}.beta, ''uni'', 0)', ...
%                       'arrayfun(@(x) x.ret_ctrl_fix, beta_all{sbj_i}.beta, ''uni'', 0)'  };
% 
% con_name_list = {'enc1-0','enc0-1','enc1','enc0', ...
%                  'ret1-0','ret0-1','ret1','ret0'};
% 
% con_eval_str_list = { 'arrayfun(@(x) x.enc1, beta_all{sbj_i}.beta, ''uni'', 0)', ...
%                       'arrayfun(@(x) x.enc0, beta_all{sbj_i}.beta, ''uni'', 0)', ...
%                       'arrayfun(@(x) x.enc1, beta_all{sbj_i}.beta_raw, ''uni'', 0)', ...
%                       'arrayfun(@(x) x.enc0, beta_all{sbj_i}.beta_raw, ''uni'', 0)', ...
%                       'arrayfun(@(x) x.ret1, beta_all{sbj_i}.beta, ''uni'', 0)', ...
%                       'arrayfun(@(x) x.ret0, beta_all{sbj_i}.beta, ''uni'', 0)', ...
%                       'arrayfun(@(x) x.ret1, beta_all{sbj_i}.beta_raw, ''uni'', 0)', ...
%                       'arrayfun(@(x) x.ret0, beta_all{sbj_i}.beta_raw, ''uni'', 0)'};

con_name_list = {'enc1-3','enc3-1','enc0-2','enc2-0', ...
                 'ret1-3','ret3-1','ret0-2','ret2-0'};

con_eval_str_list = { 'arrayfun(@(x) x.enc13, beta_all{sbj_i}.beta, ''uni'', 0)', ...
                      'arrayfun(@(x) x.enc31, beta_all{sbj_i}.beta, ''uni'', 0)', ...
                      'arrayfun(@(x) x.enc02, beta_all{sbj_i}.beta, ''uni'', 0)', ...
                      'arrayfun(@(x) x.enc20, beta_all{sbj_i}.beta, ''uni'', 0)', ...
                      'arrayfun(@(x) x.ret13, beta_all{sbj_i}.beta, ''uni'', 0)', ...
                      'arrayfun(@(x) x.ret31, beta_all{sbj_i}.beta, ''uni'', 0)', ...
                      'arrayfun(@(x) x.ret02, beta_all{sbj_i}.beta, ''uni'', 0)', ...
                      'arrayfun(@(x) x.ret20, beta_all{sbj_i}.beta, ''uni'', 0)'};

tag_nan_data = {on_on_tag|on_off_tag, on_on_tag|on_off_tag, off_on_tag|off_off_tag, off_on_tag|off_off_tag, on_on_tag|on_off_tag, on_on_tag|on_off_tag, off_on_tag|off_off_tag, off_on_tag|off_off_tag}


if opt_img == 1
    dir_glm = 'glm_norm_smooth'; is_masking = true;
elseif opt_img == 2
    dir_glm = 'glm_norm_rough'; is_masking = true;
elseif opt_img == 3
    dir_glm = 'glm_indiv_rough'; is_masking = false;
end


%
if opt_fieldmap == 1
    dir_result = 'results_spatial_202411';
    dir_fmri = sprintf('../data/data_fmri_sess%d_preprocessed',sess_i);
elseif opt_fieldmap == 2
    dir_result = 'results_spatial_unwrapped';
    dir_fmri = sprintf('../data/data_fmri_sess%d_preprocessed_corrected',sess_i);
end

%
if opt_fieldmap == 1
    if opt_img == 1
        target_reg_exp = '.*swraf.*';
    elseif opt_img == 2
        target_reg_exp = '.*\\wraf.*';
    elseif opt_img == 3
        target_reg_exp = '.*\\raf.*';
    end
elseif opt_fieldmap == 2
    if opt_img == 1
        target_reg_exp = '.*swuaf.*';
    elseif opt_img == 2
        target_reg_exp = '.*\\wuaf.*';
    elseif opt_img == 3
        target_reg_exp = '.*\\uaf.*';
    end
end


%
dir_first = sprintf('../%s/%s/%s/1st_level/_organized',dir_result, reg_setting_name, dir_glm);


flag_name_list = {'all'};
dir_second_list = {'glm_act_all'};
dir_second_list = cellfun(@(x) sprintf('../%s/%s/%s/%s', dir_result, reg_setting_name, dir_glm, x), dir_second_list, 'uni', 0);

dir_name_collect = '_organized';

% load qc data
load('../data/data_fmri_quality_control/spatial_valid_sbj.mat') % spatial_valid_sbj
% spatial_valid_sbj(17) = false;  % timing error
% save('../data/data_fmri_quality_control/spatial_valid_sbj.mat','spatial_valid_sbj')

% spatial_valid_sbj


cd(dir_working);


%% load beta - it can be more efficient (pre-loading & run), but didn't do that

beta_all = {};
parfor sbj_i = 1:length(sbj_list)
    if spatial_valid_sbj(sbj_i) == false
        beta_all{sbj_i} = [];
        continue
    end

    cd(dir_working);
    data = load(fullfile(dir_first,[sbj_list{sbj_i} '.mat']));
    beta_all{sbj_i} = data.glm;
end

%% set flag
flag_list = {};
for sbj_i = 1:length(beta_all)
    if spatial_valid_sbj(sbj_i) == false
        flag_list{sbj_i} = [];
    else
        flag_list{sbj_i} = true(1,3);
    end
end

%% save image
flag_i = 1;

for con_i = 1:length(con_name_list)
    con_name = con_name_list{con_i};

    for sbj_i = 1:length(beta_all)
        if spatial_valid_sbj(sbj_i) == false
            continue
        end

        beta = eval(con_eval_str_list{con_i});
        beta = beta(flag_list{sbj_i});
        hdr = beta_all{sbj_i}.hdr(1);
        beta = nanmean(cat(4,beta{:}),4);

        dir_temp = fullfile(dir_second_list{flag_i},dir_name_collect, con_name);
        if ~exist(dir_temp,'dir')
            mkdir(dir_temp)
        end

        % save
        niftiwrite(beta, fullfile(dir_temp, sprintf('%d.nii', sbj_i)), hdr);

    end
end


%% run second level
setting_list = combvec(1:length(con_name_list),1:length(flag_name_list));

for setting_i = 1:size(setting_list,2)
    con_i = setting_list(1,setting_i);
    flag_i = setting_list(2,setting_i);
    fprintf('\n\n\n RUNNING: con %d flag %d\n\n\n', con_i, flag_i);

    cd(dir_working);

    dir_organized = fullfile(dir_second_list{flag_i},dir_name_collect, con_name_list{con_i});
% 
%     % overall
%     flag = true(1,length(group));
%     dir_out = fullfile(dir_second_list{flag_i}, 'overall', con_name_list{con_i});
%     if ~exist(dir_out,'dir')
%         mkdir(dir_out)
%     end
%     file_list = arrayfun(@(x) fullfile(dir_organized, sprintf('%d.nii', x)), 1:length(sbj_list), 'uni', 0);
%     valid_flag = cellfun(@(x) exist(x,'file'), file_list) ~=0;
%     file_list(~valid_flag | ~flag) = [];
%     batch = jh_fmri_batch_ttest1(file_list, dir_out);
%     jh_fmri_run_batch(batch);
%     cd(dir_working);
% 
%     jh_fmri_cluster_correction(dir_out, 0.005, 10, 'pos', true);    cd(dir_working)
%     jh_fmri_cluster_correction(dir_out, 0.001, 10, 'pos', true);    cd(dir_working)


    % exp group
    flag = group == 1;
    dir_out = fullfile(dir_second_list{flag_i}, 'exp', con_name_list{con_i});
    if ~exist(dir_out,'dir')
        mkdir(dir_out)
    end
    file_list = arrayfun(@(x) fullfile(dir_organized, sprintf('%d.nii', x)), 1:length(sbj_list), 'uni', 0);
    valid_flag = cellfun(@(x) exist(x,'file'), file_list) ~=0;
    file_list(~valid_flag | ~flag) = [];
    % batch = jh_fmri_batch_ttest1(file_list, dir_out);
% 
    curr_tag = tag_nan_data{con_i};
    curr_tag((~valid_flag | ~flag)) = [];
    batch = jh_fmri_batch_ttest1(file_list(~curr_tag), dir_out);

    jh_fmri_run_batch(batch);
    cd(dir_working);

    jh_fmri_cluster_correction(dir_out, 0.005, 10, 'pos', true);    cd(dir_working)
    jh_fmri_cluster_correction(dir_out, 0.001, 10, 'pos', true);    cd(dir_working)

%     % exp group 2
%     flag = group ~= 0;
%     dir_out = fullfile(dir_second_list{flag_i}, 'exp2', con_name_list{con_i});
%     if ~exist(dir_out,'dir')
%         mkdir(dir_out)
%     end
%     file_list = arrayfun(@(x) fullfile(dir_organized, sprintf('%d.nii', x)), 1:length(sbj_list), 'uni', 0);
%     valid_flag = cellfun(@(x) exist(x,'file'), file_list) ~=0;
%     file_list(~valid_flag | ~flag) = [];
%     batch = jh_fmri_batch_ttest1(file_list, dir_out);
%     jh_fmri_run_batch(batch);
%     cd(dir_working);
% 
%     jh_fmri_cluster_correction(dir_out, 0.005, 10, 'pos', true);    cd(dir_working)
%     jh_fmri_cluster_correction(dir_out, 0.001, 10, 'pos', true);    cd(dir_working)

%     % ctrl group
%     flag = group == 0;
%     dir_out = fullfile(dir_second_list{flag_i}, 'ctrl', con_name_list{con_i});
%     if ~exist(dir_out,'dir')
%         mkdir(dir_out)
%     end
%     file_list = arrayfun(@(x) fullfile(dir_organized, sprintf('%d.nii', x)), 1:length(sbj_list), 'uni', 0);
%     valid_flag = cellfun(@(x) exist(x,'file'), file_list) ~=0;
%     file_list(~valid_flag | ~flag) = [];
%     batch = jh_fmri_batch_ttest1(file_list, dir_out);
%     jh_fmri_run_batch(batch);
%     cd(dir_working);
% 
%     jh_fmri_cluster_correction(dir_out, 0.005, 10, 'pos', true);    cd(dir_working)
%     jh_fmri_cluster_correction(dir_out, 0.001, 10, 'pos', true);    cd(dir_working)
% 
%     % exp vs. ctrl
%     dir_out = fullfile(dir_second_list{flag_i}, 'exp_vs_ctrl', con_name_list{con_i});
%     if ~exist(dir_out,'dir')
%         mkdir(dir_out)
%     end
%     file_list = arrayfun(@(x) fullfile(dir_organized, sprintf('%d.nii', x)), 1:length(sbj_list), 'uni', 0);
%     valid_flag = cellfun(@(x) exist(x,'file'), file_list) ~=0;
%     file_list1 = file_list(valid_flag & group == 1);
%     file_list2 = file_list(valid_flag & group == 0);
%     batch = jh_fmri_batch_ttest2(file_list1, file_list2, dir_out);
%     jh_fmri_run_batch(batch);
%     cd(dir_working);
% 
%     jh_fmri_cluster_correction(dir_out, 0.005, 10, 'pos', true);    cd(dir_working)
%     jh_fmri_cluster_correction(dir_out, 0.001, 10, 'pos', true);    cd(dir_working)

end


end

%%


