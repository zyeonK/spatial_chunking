%% get masks

%
sess_i = 2;
file_path_fmri = sprintf('../data/data_fmri_organized_overall_sess%d',sess_i);

fmri_info_all = {};
parfor (sbj_i = 1:length(sbj_list),7)
% for sbj_i = 1:length(sbj_list)
    data = load(fullfile(file_path_fmri,[sbj_list{sbj_i} '.mat']));
    fmri_info_all{sbj_i} = data.sbj;
end
fmri_info_all2 = fmri_info_all;

%% collect and save 

dir_working = pwd();
sbj_list = split(num2str(1:33));

TR = 2;

load('../data/data_fmri_quality_control/spatial_valid_sbj.mat') % spatial_valid_sbj

dir_glm = 'glm_indiv_rough';
dir_result = 'results_spatial_202411';  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reg_setting_name = 'task_period_basic_add_keypress'; con_type = 1;
% reg_setting_name = 'node_coin_binary'; con_type = 2;
% reg_setting_name = 'node_coin_binary_exclude_task'; con_type = 3;
reg_setting_name = 'node_boundary_2x2_target'; con_type = 4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


dir_first = sprintf('../%s/%s/%s/1st_level/_organized',dir_result, reg_setting_name, dir_glm);


%
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

[metric_all,con_name_list, roi_name_list] = func_get_beta_hpc(beta_all, fmri_info_all2, con_type);

save(sprintf('DATA_ORGANIZED_%s.mat', reg_setting_name), ...
     'metric_all','con_name_list','roi_name_list')

%% collect and save  2

dir_working = pwd();
sbj_list = split(num2str(1:33));

TR = 2;

load('../data/data_fmri_quality_control/spatial_valid_sbj.mat') % spatial_valid_sbj


dir_glm = 'glm_indiv_rough';
dir_result = 'results_spatial_202411';  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
reg_setting_name = 'task_period_basic_add_keypress'; con_type = 1;
% reg_setting_name = 'node_coin_binary'; con_type = 2;
reg_setting_name = 'node_coin_binary_exclude_task'; con_type = 3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


dir_first = sprintf('../%s/%s/%s/1st_level/_organized',dir_result, reg_setting_name, dir_glm);


%
beta_all = {};
% parfor sbj_i = 1:length(sbj_list)
for sbj_i = 1:length(sbj_list)
    if spatial_valid_sbj(sbj_i) == false
        beta_all{sbj_i} = [];
        continue
    end

    cd(dir_working);
    data = load(fullfile(dir_first,[sbj_list{sbj_i} '.mat']));
    beta_all{sbj_i} = data.glm;
end

[metric_all,con_name_list, roi_name_list] = func_get_beta_ctx(beta_all, fmri_info_all2, con_type);

save(sprintf('DATA_ORGANIZED_ctx_%s.mat', reg_setting_name), ...
     'metric_all','con_name_list','roi_name_list')

%% collect and save  3

dir_working = pwd();
sbj_list = split(num2str(1:33));

TR = 2;

load('../data/data_fmri_quality_control/spatial_valid_sbj.mat') % spatial_valid_sbj


dir_glm = 'glm_indiv_rough';
dir_result = 'results_spatial_202411';  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
reg_setting_name = 'task_period_basic_add_keypress'; con_type = 1;
% reg_setting_name = 'node_coin_binary'; con_type = 2;
% reg_setting_name = 'node_coin_binary_exclude_task'; con_type = 3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


dir_first = sprintf('../%s/%s/%s/1st_level/_organized',dir_result, reg_setting_name, dir_glm);


%
beta_all = {};
% parfor sbj_i = 1:length(sbj_list)
for sbj_i = 1:length(sbj_list)
    if spatial_valid_sbj(sbj_i) == false
        beta_all{sbj_i} = [];
        continue
    end

    cd(dir_working);
    data = load(fullfile(dir_first,[sbj_list{sbj_i} '.mat']));
    beta_all{sbj_i} = data.glm;
end

[metric_all,con_name_list, roi_name_list] = func_get_beta_hpc_ap(beta_all, fmri_info_all2, con_type);

save(sprintf('DATA_ORGANIZED_hpc_ap_%s.mat', reg_setting_name), ...
     'metric_all','con_name_list','roi_name_list')


%%
dir_working = pwd();

sbj_list = split(num2str(1:33));

TR = 2;

load('../data/data_fmri_quality_control/spatial_valid_sbj.mat') % spatial_valid_sbj




dir_result = 'results_spatial';  
dir_glm = 'glm_norm_rough';

reg_setting_name = 'task_period_basic_add_keypress';


dir_first = sprintf('../%s/%s/%s/1st_level/_organized',dir_result, reg_setting_name, dir_glm);


%
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

[metric_all, con_name_list, roi_name_list] = func_get_beta_basic(beta_all, fmri_info_all2);

save(sprintf('DATA_ORGANIZED_NORM_%s.mat', reg_setting_name), ...
     'metric_all','con_name_list','roi_name_list')

