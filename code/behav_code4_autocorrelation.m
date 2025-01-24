rng(918)

%% Set directory
data_dir = '../NC_data';
data_name  = 'group_pattern_data.mat';
%% Load data
load(fullfile(data_dir,data_name)) % avg_data, tmap_all

%% Perform autocorrelation analysis and calculate local peaks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data for autocorrelation analysis
autocorr_data = avg_data; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sample_map = ones(size(autocorr_data)); 
num_of_elements = xcorr2(sample_map);

% Filling NaN;
NaN_to_constant = 0;
autocorr_data = fillmissing(autocorr_data,'constant',NaN_to_constant);

% Compute autocorrelogram
autocorrelogram = xcorr2(autocorr_data)./num_of_elements; 

%% Visualization - autocorrelogram
figure;
hold on
imagesc(autocorrelogram');     
set(gca,'YDir','normal')
colorbar;
axis equal tight;
clim([-0.25 0.25])
set(gca,'LineWidth',1.2, 'FontWeight', 'bold','FontSize',15);
xticklabels({}); yticklabels({})
colormap('jet')
region_max_idx = imregionalmax(autocorrelogram);
[row_list, col_list] = ind2sub(size(region_max_idx), find(region_max_idx));
scatter(row_list,col_list,25,'filled','MarkerFaceColor','k')
box on
xlim([50 148])
ylim([50 148])