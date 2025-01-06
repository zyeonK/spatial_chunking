

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


%% task period
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out_path = '../data/data_regressors_spatial/task_period';
out_path2 = '../data/data_regressors_spatial/task_period_exclude_enc_ret';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist(out_path,'dir')
    mkdir(out_path)
end
if ~exist(out_path2,'dir')
    mkdir(out_path2)
end

%
for sbj_i = 1:length(sbj_list)
	sbj_name = sbj_list{sbj_i};

	out_file_name = fullfile(out_path, [sbj_name '.mat']);

	% load data
    sbj = data_all{sbj_i};
	marker = sbj.spatial.sess(end).marker;

	% initialize variable names
	start_var_list = cellfun(@(x) ['marker.' x '_start'], marker.name_all, 'UniformOutput',false);
	end_var_list = cellfun(@(x) ['marker.' x '_end'], marker.name_all, 'UniformOutput',false);

	onsets_eval_str = 'onsets = {';
	durations_eval_str = 'durations = {';
	for var_i = 1:length(marker.name_all)
	    onsets_eval_str = [onsets_eval_str start_var_list{var_i} ','];
	    durations_eval_str = [durations_eval_str end_var_list{var_i} '-' start_var_list{var_i} ','];
	end
	onsets_eval_str = [onsets_eval_str '};'];
	durations_eval_str = [durations_eval_str '};'];

	eval(onsets_eval_str);
	eval(durations_eval_str);
	names = marker.name_all;
    
    ind_order = [1 2 7 8 3 4 5 6];
    onsets = onsets(ind_order);
    durations = durations(ind_order);
    names = names(ind_order);

    % split into separate trials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    n_trial = length(onsets{1});
    n_reg = length(onsets);
    for trial_i = 2:n_trial
        onsets = [onsets, cellfun(@(x) x(trial_i), onsets(1:n_reg), 'uni', 0)];
        durations = [durations, cellfun(@(x) x(trial_i), durations(1:n_reg), 'uni', 0)];
        names = [names, names(1:n_reg)];
    end
    onsets(1:n_reg) = cellfun(@(x) x(1), onsets(1:n_reg), 'uni', 0);
    durations(1:n_reg) = cellfun(@(x) x(1), durations(1:n_reg), 'uni', 0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% save variables
	save(out_file_name, 'names','onsets','durations')
	
    % exclude encoidng and retrieval
    elim_idx = [2,6,10,14,18,22];
    onsets(elim_idx) = [];
    durations(elim_idx) = [];
    names(elim_idx) = [];
        
    save(fullfile(out_path2, [sbj_name '.mat']), 'names','onsets','durations')

    fprintf('\n%d\n',sbj_i);
end

%% nuisance regressors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out_path = '../data/data_regressors_spatial/movement';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist(out_path,'dir')
    mkdir(out_path)
end

% movement (6DoF + FD)
for sbj_i = 1:length(data_all)
    move = move_all{sbj_i}{1};  move_diff = [zeros(1,size(move,2)); diff(move)];
    rot = rot_all{sbj_i}{1};    rot_diff = [zeros(1,size(rot,2)); diff(rot)];
    fd = fd_all{sbj_i}{1};
    
    cont_reg_temp = [move, rot, fd];
    writematrix(cont_reg_temp, fullfile(out_path, sprintf('%d.txt',sbj_i)), 'delimiter','tab');
end


%% key press
% key press
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out_path1 = '../data/data_regressors_spatial/key_press_conv';
out_path2 = '../data/data_regressors_spatial/key_press_original';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist(out_path1,'dir')
    mkdir(out_path1)
end
if ~exist(out_path2,'dir')
    mkdir(out_path2)
end

for sbj_i = 1:length(data_all)
    reg_template = zeros(size(move_all{sbj_i}{1},1),1);
    
    target_list = {data_all{sbj_i}.spatial.sess(end).marker.press_left, ...
                   data_all{sbj_i}.spatial.sess(end).marker.press_right, ...
                   data_all{sbj_i}.spatial.sess(end).marker.press_front, ...
                   data_all{sbj_i}.spatial.sess(end).marker.press_confirm };
    
    regressors = {};
    regressors_original = {};
    for target_i = 1:length(target_list)
        target = target_list{target_i};
        temp = reg_template;
        for time_i = 1:length(temp)
            time_range = [(time_i-1)*TR, time_i*TR];
            if sum(target >= time_range(1) & target <= time_range(2)) ~= 0
                temp(time_i) = 1;
            end
        end
        regressors_original{target_i} = temp;
    
        temp = conv(temp, spm_hrf(TR));
        temp = temp(1:length(reg_template));
        regressors{target_i} = temp;
    end
    regressors = cell2mat(regressors);
    regressors_original = cell2mat(regressors_original);

    writematrix(regressors, fullfile(out_path1, sprintf('%d.txt',sbj_i)), 'delimiter','tab');
    writematrix(regressors_original, fullfile(out_path2, sprintf('%d.txt',sbj_i)), 'delimiter','tab');
end


%% visualize

sbj_i = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(sprintf('../data/data_regressors_spatial/task_period/%d.mat',sbj_i))
movement = readmatrix(sprintf('../data/data_regressors_spatial/movement/%d.txt',sbj_i));
keypress = readmatrix(sprintf('../data/data_regressors_spatial/key_press_conv/%d.txt',sbj_i));
keypress_orig = readmatrix(sprintf('../data/data_regressors_spatial/key_press_original/%d.txt',sbj_i));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% task period
regressors = {};
reg_template = zeros(size(movement,1),1);
for reg_i = 1:length(onsets)
    temp = reg_template;
    for onset_i = 1:length(onsets{reg_i})
        idx_start = round(onsets{reg_i}(onset_i)/TR);
        idx_end = round((onsets{reg_i}(onset_i) + durations{reg_i}(onset_i))/TR);
        temp(idx_start:idx_end) = 1;
    end
    regressors{reg_i} = temp;
end
regressors = cell2mat(regressors);

temp = conv(regressors(:,6), spm_hrf(TR));
temp = temp(1:length(regressors));
data = [regressors(:,6), temp ];

figure('position',[50 50 1800 420]);
plot(data,'linewidth',1.2)
set(gca,'LineWidth', .8,'FontSize',11, 'FontWeight','bold');
box off
ylim([-0.5 1.5])

%%% movement
% data = movement;
% figure;
% plot(data,'linewidth',1.2)
% set(gca,'LineWidth', .8,'FontSize',11, 'FontWeight','bold');
% box off


%%% key press
data = keypress;
figure('position',[50 50 1800 600]);
subplot(211)
plot(data,'linewidth',1.2)
set(gca,'LineWidth', .8,'FontSize',11, 'FontWeight','bold');
box off
subplot(212)
data = keypress_orig;
plot(data,'linewidth',1.2)
set(gca,'LineWidth', .8,'FontSize',11, 'FontWeight','bold');
box off















