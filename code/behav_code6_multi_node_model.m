rng(918)

%% Set directory
data_dir = '../data';
data_name  = 'sample_data.mat';
map_name  = 'sbj_pattern_all.mat';

%% Load data
load(fullfile(data_dir,data_name)) % table_all
load(fullfile(data_dir,map_name)) % sbj_resp_all
addpath('visualization_toolbox')

sbj_list = unique(table_all.sbj)';
num_of_sbj = length(sbj_list);


%% Multi-node modeling varying both scale and size
arena_size = 1;
step_size = 0.01;

x_coord_list=step_size:step_size:arena_size-step_size;  
y_coord_list=step_size:step_size:arena_size-step_size;   
[X,Y]=meshgrid(x_coord_list,y_coord_list);
mesh_size = size(X);

scale_list = linspace(0.15,0.5,36*arena_size);
sig_list = linspace(0.05/arena_size,0.2/arena_size,16);  

% 0-node data (random distribution without bias)
rand_data_all = {};
for data_i = 1:length(scale_list)*length(sig_list)
    rand_data = rand(mesh_size);
    rand_data_all{data_i} = (rand_data)./sum(rand_data(:));
end


corr_type = 'Pearson';
off_node_val = 0;
rho_list_0_all = {}; rho_list_4_all = {}; rho_list_5_all = {}; rho_list_9_all = {}; rho_list_16_all = {};

mdl_i = 1;
for scale_i = 1:length(scale_list)

    scale_factor = scale_list(scale_i);

    node4_center = [scale_factor/2, 1-scale_factor/2; 1-scale_factor/2, 1-scale_factor/2; scale_factor/2, scale_factor/2; 1-scale_factor/2, scale_factor/2];
    node5_center = [0.5, 0.5; scale_factor/2, 1-scale_factor/2; 1-scale_factor/2, 1-scale_factor/2; scale_factor/2, scale_factor/2; 1-scale_factor/2, scale_factor/2];
    node9_center = [0.5, 0.5; scale_factor/2, 1-scale_factor/2; 1-scale_factor/2, 1-scale_factor/2; scale_factor/2, scale_factor/2; 1-scale_factor/2, scale_factor/2; 1/2, 1-scale_factor/2; 1-scale_factor/2,1/2; 1/2, scale_factor/2; scale_factor/2, 1/2];
    node16_center = [scale_factor/2, 1-scale_factor/2; 1-scale_factor/2, 1-scale_factor/2; scale_factor/2, scale_factor/2; 1-scale_factor/2, scale_factor/2;
                   2/3-scale_factor/6,2/3-scale_factor/6; 2/3-scale_factor/6,1/3+scale_factor/6; 1/3+scale_factor/6,1/3+scale_factor/6; 1/3+scale_factor/6, 2/3-scale_factor/6;
                   1/3+scale_factor/6, 1-scale_factor/2; 2/3-scale_factor/6, 1-scale_factor/2; 1/3+scale_factor/6, scale_factor/2; 2/3-scale_factor/6,scale_factor/2; 
                   scale_factor/2, 1/3+scale_factor/6; scale_factor/2, 2/3-scale_factor/6; 1-scale_factor/2, 2/3-scale_factor/6; 1-scale_factor/2, 1/3+scale_factor/6];

    for size_i = 1:length(sig_list)
        sig = sig_list(size_i);

        %%% 0-node
        rand_data = rand_data_all{mdl_i}; 
        mdl_i = mdl_i + 1;
        rho_list = cellfun(@(x) corr(x(:), rand_data(:),'type', corr_type,'rows','pairwise'), sbj_resp_all);
        rho_list_0_all{scale_i,size_i} = [atanh(rho_list)];

        %%% 4-node
        node_center = node4_center; 
        trunc_gaussian_func = zeros(mesh_size);
        for center_i = 1:size(node_center,1)
            curr_kern = exp(-1/(2*sig^2)*((X-node_center(center_i,1)).^2+(Y-node_center(center_i,2)).^2));
            tag_zero = sqrt((X-node_center(center_i,1)).^2+(Y-node_center(center_i,2)).^2)>=sig;
            trunc_kern = curr_kern; 
            trunc_kern(tag_zero) = off_node_val;
            trunc_gaussian_func=trunc_gaussian_func+ trunc_kern;
        end
        trunc_gaussian_all = (trunc_gaussian_func)./sum(trunc_gaussian_func(:));
   
        % Pearson R
        rho_list = cellfun(@(x) corr(x(:), trunc_gaussian_all(:),'type', corr_type,'rows','pairwise'), sbj_resp_all);
        if min(pdist(node_center)) <= (2*sig)
            rho_list = NaN(1,num_of_sbj);
        end
        rho_list_4_all{scale_i,size_i} = [atanh(rho_list)];
    

        %%% 5-node
        node_center = node5_center; 
        trunc_gaussian_func = zeros(mesh_size);
        for center_i = 1:size(node_center,1)
            curr_kern = exp(-1/(2*sig^2)*((X-node_center(center_i,1)).^2+(Y-node_center(center_i,2)).^2));
            tag_zero = sqrt((X-node_center(center_i,1)).^2+(Y-node_center(center_i,2)).^2)>=sig;
            trunc_kern = curr_kern; 
            trunc_kern(tag_zero) = off_node_val;
            trunc_gaussian_func=trunc_gaussian_func+ trunc_kern;
        end
        trunc_gaussian_all = (trunc_gaussian_func)./sum(trunc_gaussian_func(:));
   
        % Pearson R
        rho_list = cellfun(@(x) corr(x(:), trunc_gaussian_all(:),'type', corr_type,'rows','pairwise'), sbj_resp_all);
        if min(pdist(node_center)) <= (2*sig)
            rho_list = NaN(1,num_of_sbj);
        end
        rho_list_5_all{scale_i,size_i} = [atanh(rho_list)];
    

        %%% 9-node
        node_center = node9_center; 
        trunc_gaussian_func = zeros(mesh_size);
        for center_i = 1:size(node_center,1)
            curr_kern = exp(-1/(2*sig^2)*((X-node_center(center_i,1)).^2+(Y-node_center(center_i,2)).^2));
            tag_zero = sqrt((X-node_center(center_i,1)).^2+(Y-node_center(center_i,2)).^2)>=sig;
            trunc_kern = curr_kern; 
            trunc_kern(tag_zero) = off_node_val;
            trunc_gaussian_func=trunc_gaussian_func+ trunc_kern;
        end
        trunc_gaussian_all = (trunc_gaussian_func)./sum(trunc_gaussian_func(:));
   
        % Pearson R
        rho_list = cellfun(@(x) corr(x(:), trunc_gaussian_all(:),'type', corr_type,'rows','pairwise'), sbj_resp_all);
        if min(pdist(node_center)) <= (2*sig)
            rho_list = NaN(1,num_of_sbj);
        end
        rho_list_9_all{scale_i,size_i} = [atanh(rho_list)];
    

        %%% 16-node
        node_center = node16_center; 
        trunc_gaussian_func = zeros(mesh_size);
        for center_i = 1:size(node_center,1)
            curr_kern = exp(-1/(2*sig^2)*((X-node_center(center_i,1)).^2+(Y-node_center(center_i,2)).^2));
            tag_zero = sqrt((X-node_center(center_i,1)).^2+(Y-node_center(center_i,2)).^2)>=sig;
            trunc_kern = curr_kern; 
            trunc_kern(tag_zero) = off_node_val;
            trunc_gaussian_func=trunc_gaussian_func+ trunc_kern;
        end
        trunc_gaussian_all = (trunc_gaussian_func)./sum(trunc_gaussian_func(:));
   
        % Pearson R
        rho_list = cellfun(@(x) corr(x(:), trunc_gaussian_all(:),'type', corr_type,'rows','pairwise'), sbj_resp_all);
        if min(pdist(node_center)) <= (2*sig)
            rho_list = NaN(1,num_of_sbj);
        end
        rho_list_16_all{scale_i,size_i} = [atanh(rho_list)];
    
        
    end
end

% save data
save([data_dir,'multi_node_model_fit.mat'], 'rho_list_0_all', 'rho_list_4_all','rho_list_5_all','rho_list_9_all','rho_list_16_all')

%% Find best fitting model based on correlation coefficient
% spacing scale and node sizes

max_idx = []; max_val = []; 
node_n_num = 5; 

metric_all = {rho_list_0_all, rho_list_4_all, rho_list_5_all, rho_list_9_all, rho_list_16_all};
for node_i = 1:node_n_num
    curr_node = metric_all{node_i};
    for sbj_i  = 1:num_of_sbj
        curr_sbj_data = cell2mat(cellfun(@(x) ~isnan(x(sbj_i))*x(sbj_i), curr_node,'UniformOutput',0));
        [max_val(node_i, sbj_i),I] = max(curr_sbj_data(:));
        [I_scale, I_size] = ind2sub([length(scale_list), length(sig_list)],I);
        max_idx{node_i}(sbj_i,:) = [I_scale, I_size];
    end

end

max_rho_all = max_val; % model fitness (correlation coeff.)

%% Visualization
% Histogram: number of participants whose responses most highly correlated with each chunking model

[M,I] = max(max_rho_all);
figure;
hold on
X = [sum(I == 1),sum(I == 2),sum(I == 3),sum(I == 4),sum(I == 5)];
set(gca,'LineWidth',1.5, 'FontWeight', 'bold','FontSize',15);
jh_bar(X)
xticks([1:5])
xticklabels({'0', '4', '5', '9','16'})
ylim([0 num_of_sbj])

%% Statistical tests

% multi-node-modeling results (model fitness)
p=[];
[h,p(1),~,stats] = ttest(max_rho_all(1,:));
fprintf('t-test 0-node: p = %.3f, t(%d) = %.3f\n', p(1), stats.df, stats.tstat)
[h,p(2),~,stats] = ttest(max_rho_all(2,:));
fprintf('t-test 4-node: p = %.3f, t(%d) = %.3f\n', p(2), stats.df, stats.tstat)
[h,p(3),~,stats] = ttest(max_rho_all(3,:));
fprintf('t-test 5-node: p = %.3f, t(%d) = %.3f\n', p(3), stats.df, stats.tstat)
[h,p(4),~,stats] = ttest(max_rho_all(4,:));
fprintf('t-test 9-node: p = %.3f, t(%d) = %.3f\n', p(4), stats.df, stats.tstat)
[h,p(5),~,stats] = ttest(max_rho_all(5,:));
fprintf('t-test 16-node: p = %.3f, t(%d) = %.3f\n', p(5), stats.df, stats.tstat)
[a,b,c,d] = fdr_bh(p)

% multi-node-modeling results (comparison with 0-node model)
p=[];
[h,p(1),~,stats] = ttest(max_rho_all(1,:),max_rho_all(2,:));
fprintf('0-node t-test with 4-node: p = %.3f, t(%d) = %.3f\n', p(1), stats.df, stats.tstat)
[h,p(2),~,stats] = ttest(max_rho_all(1,:),max_rho_all(3,:));
fprintf('0-node t-test with 5-node: p = %.3f, t(%d) = %.3f\n', p(2), stats.df, stats.tstat)
[h,p(3),~,stats] = ttest(max_rho_all(1,:),max_rho_all(4,:));
fprintf('0-node t-test with 9-node: p = %.3f, t(%d) = %.3f\n', p(3), stats.df, stats.tstat)
[h,p(4),~,stats] = ttest(max_rho_all(1,:),max_rho_all(5,:));
fprintf('0-node t-test with 16-node: p = %.3f, t(%d) = %.3f\n', p(4), stats.df, stats.tstat)
[a,b,c,d] = fdr_bh(p)

% multi-node modeling results (comparison with 9-node model)
nine_node_idx = 4;
p=[];
[h,p(1),~,stats] = ttest(max_rho_all(nine_node_idx,:),max_rho_all(1,:));
fprintf('9-node t-test with 0-node: p = %.3f, t(%d) = %.3f\n', p(1), stats.df, stats.tstat)
[h,p(2),~,stats] = ttest(max_rho_all(nine_node_idx,:),max_rho_all(2,:));
fprintf('9-node t-test with 4-node: p = %.3f, t(%d) = %.3f\n', p(1), stats.df, stats.tstat)
[h,p(3),~,stats] = ttest(max_rho_all(nine_node_idx,:),max_rho_all(3,:));
fprintf('9-node t-test with 5-node: p = %.3f, t(%d) = %.3f\n', p(2), stats.df, stats.tstat)
[h,p(4),~,stats] = ttest(max_rho_all(nine_node_idx,:),max_rho_all(5,:));
fprintf('9-node t-test with 16-node: p = %.3f, t(%d) = %.3f\n', p(3), stats.df, stats.tstat)
[a,b,c,d] = fdr_bh(p)

%% Define each individual's best fitting 9-node model
node_i = 4; % 9-node model's index

x_coord_list=step_size:step_size:arena_size-step_size;  
y_coord_list=step_size:step_size:arena_size-step_size;   
[X,Y]=meshgrid(x_coord_list,y_coord_list);

sbj_node_model = {}; node_center = {}; node_sizes = {};
for sbj_i = 1:num_of_sbj

    scale_factor = scale_list(max_idx{node_i}(sbj_i,1)); 
    sig = sig_list(max_idx{node_i}(sbj_i,2));

    node_center = [0.5, 0.5; scale_factor/2, 1-scale_factor/2; 1-scale_factor/2, 1-scale_factor/2; scale_factor/2, scale_factor/2; 1-scale_factor/2, scale_factor/2; 1/2, 1-scale_factor/2; 1-scale_factor/2,1/2; 1/2, scale_factor/2; scale_factor/2, 1/2];

    curr_sbj_model = zeros(mesh_size);
    for center_i = 1:size(node_center,1)
        curr_kern = exp(-1/(2*sig^2)*((X-node_center(center_i,1)).^2+(Y-node_center(center_i,2)).^2));
        tag_zero = sqrt((X-node_center(center_i,1)).^2+(Y-node_center(center_i,2)).^2)>=sig;
        trunc_kern = curr_kern; trunc_kern(tag_zero) = 0;
        curr_sbj_model =curr_sbj_model+ trunc_kern;
    end

    sbj_node_model{sbj_i} = curr_sbj_model;
    node_centers{sbj_i} = node_center;
    node_sizes{sbj_i} = sig;
end

save([data_dir,'/multi_node_model_info.mat'], 'node_centers', 'node_sizes')
