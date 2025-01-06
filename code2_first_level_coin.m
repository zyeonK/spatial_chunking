%% parameters

opt_img = 1; % norm_smooth / norm_rough / indiv_rough
opt_fieldmap = 1; % no-correction / correction

sess_i = 2;

opt_list = { [1 1], [2 1], [3 1] };



for opt_i = 1:length(opt_list)
    opt_img = opt_list{opt_i}(1);
    opt_fieldmap = opt_list{opt_i}(2);



%% settings
sbj_list = split(num2str(1:33));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reg_setting_name = 'task_period_basic_add_keypress';
% reg_setting_name = 'task_period_basic';
% reg_setting_name = 'node_coin_binary';
% reg_setting_name = 'node_boundary_continuous_each';
% reg_setting_name = 'boundary_node_continuous_target';
reg_setting_name = 'node_onoff_boundary_dist_each';

% dir_reg_binary =  '../results_spatial_202411/data_regressors_spatial/node_coin_binary_conv';
% dir_reg_binary =  '../results_spatial_202411/data_regressors_spatial/boundary_node_continuous_target_conv';
dir_reg_binary =  '../results_spatial_202411/data_regressors_spatial/node_onoff_boundary_dist_each_conv';
dir_reg =  '../data/data_regressors_spatial/task_period';
dir_reg_movement = '../data/data_regressors_spatial/movement';
dir_reg_keypress = '../data/data_regressors_spatial/key_press_conv';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
dir_working = pwd();
path_behav = '../data/data_processed_sbj';

%
% con_names = {'enc','ret','ctrl_enc','ctrl_ret','enc-ctrl','ret-ctrl', 'enc-ctrl2', 'ret-ctrl2'};
con_names = {'enc','ret'};  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% intentionally skipped here

%
dir_list_func = {'spatial'};
 

%
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
dir_first = sprintf('../%s/%s/%s/1st_level',dir_result, reg_setting_name, dir_glm);
dir_organized = sprintf('../%s/%s/%s/1st_level/_organized', dir_result, reg_setting_name, dir_glm);


% load qc data
load('../data/data_fmri_quality_control/spatial_valid_sbj.mat') % spatial_valid_sbj
% spatial_valid_sbj(17) = false;  % timing error
% save('../data/data_fmri_quality_control/spatial_valid_sbj.mat','spatial_valid_sbj') % spatial_valid_sbj


%% run first-level analysis

err_list = zeros(1,length(sbj_list));
parfor sbj_i = 1:length(sbj_list)
    if spatial_valid_sbj(sbj_i) == false
        continue
    end

    try
        cd(dir_working);
    
        sbj_name = sbj_list{sbj_i};
    
        % output directory
        dir_out = fullfile(dir_first, sbj_name);
        dir_reg_out = fullfile(dir_first, sbj_name, 'reg');
        if ~exist(dir_out,'dir')
            mkdir(dir_out)
        end
        if ~exist(dir_reg_out, 'dir')
            mkdir(dir_reg_out)
        end
    
        % get functional & main regressor
        in_dir = fullfile(dir_fmri, sbj_name, dir_list_func{1});
        func = {jh_get_file_list(in_dir, target_reg_exp,'regexp')};
        task_reg = {fullfile(dir_reg, sprintf('%s.mat',sbj_name))};

    
        % get regressors and save
        data0 = readmatrix(fullfile(dir_reg_binary, sprintf('%s.txt',sbj_name)));
        data1 = readmatrix(sprintf('%s/%d.txt',dir_reg_movement,sbj_i)); temp = data1;
        data2 = readmatrix(sprintf('%s/%d.txt',dir_reg_keypress,sbj_i));     %%%%%%%%%%%%%%%%%%
        temp = [data0, data1, data2];

        writematrix(temp, sprintf('%s/regressor.txt',dir_reg_out))
        cont_reg = {sprintf('%s/regressor.txt',dir_reg_out)};
    
        % run
        batch = jh_fmri_first_level_analysis({dir_out}, func, task_reg, cont_reg, [], [], is_masking, 'none');
        jh_fmri_run_batch(batch);
    catch
        err_list(sbj_i) = 1;
    end
        cd(dir_working);

end



%% collect data

if ~exist(dir_organized,'dir')
    mkdir(dir_organized)
end

for sbj_i = 1:length(sbj_list)
    if spatial_valid_sbj(sbj_i) == false
        continue
    end

   
    sbj_name = sbj_list{sbj_i};

    clear glm beta_raw  beta  hdr

    glm = struct();
    for trial_i = 1:3
        temp_raw = struct();
        temp = struct();

        % set file names
        file_name_enc_fix = sprintf('beta_%04d.nii', 1+(trial_i-1)*8);
        file_name_enc = sprintf('beta_%04d.nii', 2+(trial_i-1)*8);
        file_name_ret_fix = sprintf('beta_%04d.nii', 5+(trial_i-1)*8);
        file_name_ret = sprintf('beta_%04d.nii', 6+(trial_i-1)*8);
        file_name_ctrl_fix = sprintf('beta_%04d.nii', 7+(trial_i-1)*8);
        file_name_ctrl = sprintf('beta_%04d.nii', 8+(trial_i-1)*8);

        file_name_enc1 = sprintf('beta_%04d.nii', 1 + 8*3 + (trial_i-1)*2);
        file_name_enc0 = sprintf('beta_%04d.nii', 2 + 8*3 + (trial_i-1)*2);

        file_name_ret1 = sprintf('beta_%04d.nii', 1 + 8*3 + 2*3 + (trial_i-1)*2);
        file_name_ret0 = sprintf('beta_%04d.nii', 2 + 8*3 + 2*3 +(trial_i-1)*2);


        in_path = fullfile(dir_first, sbj_name);
        func_load_file = @(fname) double(niftiread(fullfile(in_path, fname)));
        temp_raw.enc = func_load_file(file_name_enc);
        temp_raw.ret = func_load_file(file_name_ret);
        temp_raw.enc_fix = func_load_file(file_name_enc_fix);
        temp_raw.ret_fix = func_load_file(file_name_ret_fix);

        temp_raw.ctrl = func_load_file(file_name_ctrl);
        temp_raw.ctrl_fix = func_load_file(file_name_ctrl_fix);

        temp_raw.enc1 = func_load_file(file_name_enc1);
        temp_raw.enc0 = func_load_file(file_name_enc0);
        temp_raw.ret1 = func_load_file(file_name_ret1);
        temp_raw.ret0 = func_load_file(file_name_ret0);

        temp_list = dir(fullfile(in_path, 'beta*'));
        temp_raw.baseline = func_load_file(sprintf('beta_%04d.nii', length(temp_list)));

        % contrast images
        temp.enc = temp_raw.enc - temp_raw.enc_fix;
        temp.ret = temp_raw.ret - temp_raw.ret_fix;

        temp.enc1 = temp_raw.enc1 - temp_raw.enc0;
        temp.enc0 = temp_raw.enc0 - temp_raw.enc1;
        temp.ret1 = temp_raw.ret1 - temp_raw.ret0;
        temp.ret0 = temp_raw.ret0 - temp_raw.ret1;

        temp.enc_ctrl = temp_raw.enc - temp_raw.ctrl;
        temp.ret_ctrl = temp_raw.ret - temp_raw.ctrl;

        temp.enc_ctrl_fix = (temp_raw.enc - temp_raw.enc_fix) - (temp_raw.ctrl - temp_raw.ctrl_fix);
        temp.ret_ctrl_fix = (temp_raw.ret - temp_raw.ret_fix) - (temp_raw.ctrl - temp_raw.ctrl_fix);

        % save
        beta_raw(trial_i) = temp_raw;
        beta(trial_i) = temp;

        temp = niftiinfo(fullfile(in_path, file_name_enc_fix));
        temp.Datatype = 'double';
        temp.Description = '';
        hdr(trial_i) = temp;
    end

    
    glm.beta_raw = beta_raw;
    glm.beta = beta;
    glm.hdr = hdr;
    
    save(sprintf('%s/%s.mat',dir_organized, sbj_list{sbj_i}), 'glm');
    fprintf('\n%d\n',sbj_i)    
end


%%

%%




end







































