%% parameters

opt_img = 1; % norm_smooth / norm_rough / indiv_rough
opt_fieldmap = 1; % no-correction / correction

sess_i = 2;

opt_list = { [1 1] }; % {[option: opt_img, opt_fieldmap]}

for opt_i = 1:length(opt_list)
    opt_img = opt_list{opt_i}(1);
    opt_fieldmap = opt_list{opt_i}(2);


%% settings
sbj_list = split(num2str(1:33));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
reg_setting_name = 'task_period_basic';
reg_setting_name = 'task_period_basic_add_keypress';
% reg_setting_name = 'prepost_coin_binary'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
dir_working = pwd();
path_behav = '../data/data_processed_sbj';


con_name_list = {'enc_ctrl','ret_ctrl'};


con_eval_str_list = { 'arrayfun(@(x) x.enc_ctrl, beta_all{sbj_i}.beta, ''uni'', 0)', ...
                      'arrayfun(@(x) x.ret_ctrl, beta_all{sbj_i}.beta, ''uni'', 0)' };


% con_name_list = {'enc_raw', 'ret_raw', 'enc_fix', 'ret_fix','enc_ctrl','ret_ctrl','enc_ctrl_fix','ret_ctrl_fix'};
% 
% 
% con_eval_str_list = { 'arrayfun(@(x) x.enc, beta_all{sbj_i}.beta_raw, ''uni'', 0)', ...
%                       'arrayfun(@(x) x.ret, beta_all{sbj_i}.beta_raw, ''uni'', 0)', ...
%                       'arrayfun(@(x) x.enc, beta_all{sbj_i}.beta, ''uni'', 0)', ...
%                       'arrayfun(@(x) x.ret, beta_all{sbj_i}.beta, ''uni'', 0)', ...
%                       'arrayfun(@(x) x.enc_ctrl, beta_all{sbj_i}.beta, ''uni'', 0)', ...
%                       'arrayfun(@(x) x.ret_ctrl, beta_all{sbj_i}.beta, ''uni'', 0)', ...
%                       'arrayfun(@(x) x.enc_ctrl_fix, beta_all{sbj_i}.beta, ''uni'', 0)', ...
%                       'arrayfun(@(x) x.ret_ctrl_fix, beta_all{sbj_i}.beta, ''uni'', 0)'  };

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
dir_first = sprintf('../%s/%s/%s/1st_level/_organized',dir_result, reg_setting_name, dir_glm);


flag_name_list = {'all'};
dir_second_list = {'glm_act_all'};
dir_second_list = cellfun(@(x) sprintf('../%s/%s/%s/%s', dir_result, reg_setting_name, dir_glm, x), dir_second_list, 'uni', 0);

dir_name_collect = '_organized';
dir_name_correlation = '_correlation';

% load qc data
load('../data/data_fmri_quality_control/spatial_valid_sbj.mat') % spatial_valid_sbj
% spatial_valid_sbj(17) = false;  % timing error
% save('../data/data_fmri_quality_control/spatial_valid_sbj.mat','spatial_valid_sbj')

% spatial_valid_sbj


cd(dir_working);

%% set behavior metrics - OLD VERSION

%%%%%%%%%%%%%%%%%%%%%%%%%% refer ../code_em1_behavior/code5_metrics_correlation

% % accuracy and training metric
% aux_func_get_residual = @(mdl) mdl.Residuals.Raw;
% func_get_residual = @(x,y) aux_func_get_residual(fitlm(x(:),y(:)));
% 
% sess1_behav_spatial_list = {spatial_metric_sess_wise{1}, spatial_metric_sess_wise_25{1}, spatial_metric_sess_wise_50{1}, spatial_metric_sess_wise_75{1} };
% sess8_behav_spatial_list = {spatial_metric_sess_wise{8}, spatial_metric_sess_wise_25{8}, spatial_metric_sess_wise_50{8}, spatial_metric_sess_wise_75{8} };
% sess8_behav_spatial_raw_list = {spatial_metric_sess_wise_raw{8}, spatial_metric_sess_wise_raw{8}, spatial_metric_sess_wise_raw{8}, spatial_metric_sess_wise_raw{8} };
% training_behav_spatial_list = {spatial_metric_training, spatial_metric_training_25, spatial_metric_training_50, spatial_metric_training_75};
% 
% training_cov_behav_spatial_list = {};
% template = nan(1,length(sess1_behav_spatial_list{1}));
% for i = 1:length(training_behav_spatial_list)
%     training_cov_behav_spatial_list{i} = template;
%     training_cov_behav_spatial_list{i}(group==1) = ...
%         func_get_residual(sess1_behav_spatial_list{i}(group==1), training_behav_spatial_list{i}(group==1));
% end
% 
% sess8_behav_spatial_raw_name_list = {'sess8_spatial_raw','sess8_spatial_raw_25','sess8_spatial_raw_50','sess8_spatial_raw_75'};
% sess8_behav_spatial_name_list = {'sess8_spatial_acc','sess8_spatial_acc_25','sess8_spatial_acc_50','sess8_spatial_acc_75'};
% training_behav_spatial_name_list = {'training_spatial','training_spatial_acc_25','training_spatial_acc_50','training_spatial_acc_75'};
% training_cov_behav_spatial_name_list = {'training_cov_spatial','training_cov_spatial_acc_25','training_cov_spatial_acc_50','training_cov_spatial_acc_75'};
% 
% 
% % chunking index
% load('chunking_index.mat')
% chunking_idx1 = nan(1,length(sbj_list));
% chunking_idx2 = chunking_idx1;
% chunking_idx1(group~=0) = chunking_all(1,:);
% chunking_idx2(group~=0) = chunking_all(2,:);
% 
% chunking_idx1_cov_raw = func_get_residual(spatial_metric_sess_wise_raw{8}, chunking_idx1);
% chunking_idx2_cov_raw = func_get_residual(spatial_metric_sess_wise_raw{8}, chunking_idx2);
% chunking_idx1_cov_acc = func_get_residual(spatial_metric_sess_wise{8}, chunking_idx1);
% chunking_idx2_cov_acc = func_get_residual(spatial_metric_sess_wise{8}, chunking_idx2);
% 
% chunking_idx_list = {chunking_idx1, chunking_idx2, chunking_idx1_cov_raw, chunking_idx2_cov_raw, ...
%                      chunking_idx1_cov_acc, chunking_idx2_cov_acc };
% chunking_index_name_list = {'chunking_idx_trunc_pixcorr','chunking_idx_trunc_snr', ...
%                             'chunking_idx_trunc_pixcorr_cov_raw','chunking_idx_trunc_snr_cov_raw', ...
%                             'chunking_idx_trunc_pixcorr_cov_acc','chunking_idx_trunc_snr_cov_acc'};
% 
% % combine
% behav_metric_list = [sess8_behav_spatial_list, sess8_behav_spatial_raw_list, ...
%                     training_behav_spatial_list, training_cov_behav_spatial_list, ...
%                     chunking_idx_list ];
% behav_metric_name_list = [ sess8_behav_spatial_name_list, sess8_behav_spatial_raw_name_list, ...
%                            training_behav_spatial_name_list, training_cov_behav_spatial_name_list, ...
%                            chunking_index_name_list ];

% behav_metric_name_list = cellfun(@(x) ['resp_',x], behav_metric_name_list, 'uni', 0);  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% set behavior metrics - OLD VERSION

%%%%%%%%%%%%%%%%%%%%%%%%%% refer ../code_em1_behavior/code5_metrics_correlation

% chunking index
load('chunking_idx.mat')

tag_behav = 4;

behav_metric_list = arrayfun(@(i) data(i,:), 1:size(data,1), uni=0);
behav_metric_name_list = arrayfun(@(i) sprintf('chunking%d',i),1:size(data,1), uni=0);

behav_metric_list = behav_metric_list(tag_behav)
behav_metric_name_list = behav_metric_name_list(tag_behav)
behav_metric_list{2} = fmri_chunking(:,1);
behav_metric_name_list{2} = 'first_half';
behav_metric_list{2} = fmri_chunking(:,2);
behav_metric_name_list{2} = 'second_half';
%% load beta - it can be more efficient (pre-loading & run), but didn't do that

beta_all = {};
for sbj_i = 1:length(sbj_list)
    if spatial_valid_sbj(sbj_i) == false
        beta_all{sbj_i} = [];
        continue
    end

    cd(dir_working);
    data = load(fullfile(dir_first,[sbj_list{sbj_i} '.mat']));
    beta_all{sbj_i} = data.glm;
end

%% cleanup          MUST RUN THIS IF YOU RE-RUN THE ANALYSES
    
% behav_flag_name_list = {'overall','exp','ctrl'};
% for con_i = 1:length(con_name_list)
%     for behav_flag_i = 1:length(behav_flag_name_list)
%            dir_out = fullfile(dir_second_list{1}, behav_flag_name_list{behav_flag_i}, con_name_list{con_i}, dir_name_correlation);
%            rmdir(dir_out,'s')
%     end
% end


%% run second level
setting_list = combvec(1:length(con_name_list),1:length(flag_name_list));

for setting_i = 1:size(setting_list,2)
    con_i = setting_list(1,setting_i);
    flag_i = setting_list(2,setting_i);
    fprintf('\n\n\n RUNNING: con %d flag %d\n\n\n', con_i, flag_i);

    cd(dir_working);

    dir_organized = fullfile(dir_second_list{flag_i},dir_name_collect, con_name_list{con_i});

    % run

%     behav_flag_list = {true(1,length(group)), group==1, group~=0, group==0};
%     behav_flag_name_list = {'overall','exp', 'exp2','ctrl'};
    behav_flag_list = {group==1};
    behav_flag_name_list = {'exp'};

    for behav_flag_i = 1:length(behav_flag_list)
        behav_flag = behav_flag_list{behav_flag_i};
        behav_flag_name = behav_flag_name_list{behav_flag_i};
    
        dir_organized = fullfile(dir_second_list{flag_i},dir_name_collect, con_name_list{con_i});
        
        file_list = arrayfun(@(x) fullfile(dir_organized, sprintf('%d.nii', x)), find(behav_flag), 'uni', 0);
        valid_flag_neural = cellfun(@(x) exist(x,'file'), file_list) ~=0;
        file_list = file_list(valid_flag_neural);
        beta_list = cellfun(@(x) niftiread(x), file_list, 'uni',0);
        
        beta = cat(4,beta_list{:});
        beta_flat = reshape(beta, prod(size(beta,1:3)), size(beta,4));
    
        for behav_i = 1:length(behav_metric_list)
            % set directory
            dir_name_behav = behav_metric_name_list{behav_i};
    
            dir_out = fullfile(dir_second_list{flag_i}, behav_flag_name, con_name_list{con_i}, dir_name_correlation, dir_name_behav);
            if ~exist(dir_out,'dir')
                mkdir(dir_out)
            end
    
            % get behavioral
%             behav = behav_metric_list{behav_i}(behav_flag);
            behav = behav_metric_list{behav_i};
            behav_name = behav_metric_name_list{behav_i};
            behav(~valid_flag_neural) = [];
            behav = behav(:);
    
            flag = ~isnan(behav);
    
            % run
            temp = fopen(fullfile(dir_out, sprintf('__n = %d',sum(flag))),'w'); fclose(temp);
            if sum(flag) >= 7
                % correlation 
                [r,p] = corr(beta_flat(:,flag)', behav(flag));
                r = reshape(r, size(beta,1:3));
                p = reshape(p, size(beta,1:3));
                hdr = niftiinfo(file_list{1});
                hdr.Datatype = 'double';

                % save
                temp_fid = fopen(fullfile(dir_working, dir_out, '_file_list.txt'),'w'); 
                temp_list = cellfun(@(x) fullfile(dir_working, x), file_list(flag), 'uni', 0); 
                for i = 1:length(temp_list); fprintf(temp_fid, '%s\n', temp_list{i}); end
                fclose(temp_fid);
                
                temp_fid = fopen(fullfile(dir_working, dir_out, '_behav.txt'),'w'); 
                temp_list = behav(flag);
                for i = 1:length(temp_list); fprintf(temp_fid, '%f\n', temp_list(i)); end
                fclose(temp_fid);
                

                niftiwrite(r, fullfile(dir_out, '_r_raw'), hdr,'compressed',true)
                niftiwrite(-log10(p), fullfile(dir_out, '_p_raw'), hdr,'compressed',true)

                temp = -log10(p); temp(r<0) = -temp(r<0); 
                niftiwrite(temp, fullfile(dir_out, '_p_raw_signed'), hdr,'compressed',true)

                temp = r; temp(p<0.05) = 0; niftiwrite(temp, fullfile(dir_out, '_r_cutoff_0_05'), hdr,'compressed',true)
                temp = r; temp(p<0.01) = 0; niftiwrite(temp, fullfile(dir_out, '_r_cutoff_0_01'), hdr,'compressed',true)
                temp = r; temp(p<0.1) = 0; niftiwrite(temp, fullfile(dir_out, '_r_cutoff_0_1'), hdr,'compressed',true)
                temp = r; temp(p>=0.05 | r<0) = 0; niftiwrite(temp, fullfile(dir_out, '_r_cutoff_0_05_pos'), hdr,'compressed',true)
                temp = r; temp(p>=0.01 | r<0) = 0; niftiwrite(temp, fullfile(dir_out, '_r_cutoff_0_01_pos'), hdr,'compressed',true)
                temp = r; temp(p>=0.1 | r<0) = 0; niftiwrite(temp, fullfile(dir_out, '_r_cutoff_0_1_pos'), hdr,'compressed',true)
                temp = r; temp(p>=0.05 | r>0) = 0; niftiwrite(temp, fullfile(dir_out, '_r_cutoff_0_05_neg'), hdr,'compressed',true)
                temp = r; temp(p>=0.01 | r>0) = 0; niftiwrite(temp, fullfile(dir_out, '_r_cutoff_0_01_neg'), hdr,'compressed',true)
                temp = r; temp(p>=0.1 | r>0) = 0; niftiwrite(temp, fullfile(dir_out, '_r_cutoff_0_1_neg'), hdr,'compressed',true)


                cd(dir_working);
                batch = jh_fmri_batch_anova1(file_list(flag), dir_out, {behav(flag)}, {behav_name});
                jh_fmri_run_batch(batch);
                cd(dir_working);
                
                data = load(fullfile(dir_out, 'SPM.mat'));
                t = double(niftiread(fullfile(dir_out, 'spmF_0001.nii')));
                hdr = niftiinfo(fullfile(dir_out, 'spmF_0001.nii')); 
                hdr.Datatype = 'double';
                mask = niftiread(fullfile(dir_out, 'mask.nii')) == 1;
                
                p = zeros(size(t));
                p(t>0 & mask) = 1 - fcdf(t(t>0),1,data.SPM.xX.erdf);
                p(t<0 & mask) = 1 - fcdf(-t(t<0),1,data.SPM.xX.erdf);
                p(t==0 | ~mask ) = nan;
                p_fwe = min(p*sum(~isnan(p(:))), 1);
                p_fdr = spm_P_FDR(p);
                
                p = -log10(p);
                p_fwe = -log10(p_fwe);
                p_fdr = -log10(p_fdr);
                
                p(t<0) = -p(t<0);
                p_fwe(t<0) = -p_fwe(t<0);
                p_fdr(t<0) = -p_fdr(t<0);
                
                niftiwrite(p,fullfile(dir_out,'spm_p'),hdr,'compressed',true)
                niftiwrite(p_fwe,fullfile(dir_out,'spm_p_fwe'),hdr,'compressed',true)
                niftiwrite(p_fdr,fullfile(dir_out,'spm_p_fdr'),hdr,'compressed',true)

                jh_fmri_cluster_correction(dir_out, 0.005, 10, 'pos', true);    cd(dir_working)
                jh_fmri_cluster_correction(dir_out, 0.001, 10, 'pos', true);    cd(dir_working)

            else
                niftiwrite(nan, fullfile(dir_out, '_r_raw'), 'compressed',true)
                niftiwrite(nan, fullfile(dir_out, '_p_raw'), 'compressed',true)
                niftiwrite(nan, fullfile(dir_out, '_p_raw_signed'), 'compressed',true)
                niftiwrite(nan, fullfile(dir_out, '_r_cutoff_0_05'), 'compressed',true)
                niftiwrite(nan, fullfile(dir_out, '_r_cutoff_0_01'), 'compressed',true)
                niftiwrite(nan, fullfile(dir_out, '_r_cutoff_0_1'), 'compressed',true)
                niftiwrite(nan, fullfile(dir_out, '_r_cutoff_0_05_pos'), 'compressed',true)
                niftiwrite(nan, fullfile(dir_out, '_r_cutoff_0_01_pos'), 'compressed',true)
                niftiwrite(nan, fullfile(dir_out, '_r_cutoff_0_1_pos'), 'compressed',true)
                niftiwrite(nan, fullfile(dir_out, '_r_cutoff_0_05_neg'), 'compressed',true)
                niftiwrite(nan, fullfile(dir_out, '_r_cutoff_0_01_neg'), 'compressed',true)
                niftiwrite(nan, fullfile(dir_out, '_r_cutoff_0_1_neg'), 'compressed',true)
            end
        end
    end
    %



end

%%




end

%%


