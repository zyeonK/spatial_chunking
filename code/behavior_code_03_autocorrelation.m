clear all
rng(918)

%% Set directory
data_dir = '../data';
data_name  = 'group_pattern_data.mat';

%% Load data
load(fullfile(data_dir,data_name)) % avg_data, tmap_all

%% Perform autocorrelation analysis and calculate local peaks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data for autocorrelation analysis
autocorr_data = avg_data; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mu = mean(autocorr_data(:),'omitnan');  sig = std(autocorr_data(:),'omitnan');
autocorr_data = (autocorr_data - mu) / sig;               % z-score
autocorr_data = fillmissing(autocorr_data,'constant',0); 

sample_map = ones(size(autocorr_data)); 
num_of_elements = xcorr2(sample_map);

% Compute autocorrelogram
autocorrelogram = xcorr2(autocorr_data)./num_of_elements; 

%% Visualization - autocorrelogram (Fig 2c)
figure;
hold on
imagesc(autocorrelogram');     
set(gca,'YDir','normal')
colorbar;
axis equal tight;
clim([-0.6 0.6])
set(gca,'LineWidth',1.2, 'FontWeight', 'normal','FontSize',15);
xticklabels({}); yticklabels({})
colormap('jet')
region_max_idx = imregionalmax(autocorrelogram');
[H,W] = size(region_max_idx);
ctr = [ceil(H/2), ceil(W/2)];
pad = 1;
region_max_idx(ctr(1)-pad:ctr(1)+pad, ctr(2)-pad:ctr(2)+pad) = false;

[row_list, col_list] = ind2sub(size(region_max_idx), find(region_max_idx));
scatter(col_list,row_list,25,'filled','MarkerFaceColor','k')
box on
xlim([50 148])
ylim([50 148])

