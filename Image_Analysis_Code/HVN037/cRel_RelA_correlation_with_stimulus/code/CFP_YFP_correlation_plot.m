% load('/Users/yijiachen/Documents/B_Cell_project/HVN037/data_saved/HVN037_CFP_YFP_nuc_median_mean_int.mat')
% load('/Users/yijiachen/Documents/B_Cell_project/HVN037/data_saved/HVN037_CFP_YFP_wc_median_mean_int.mat')
load('/Users/yijiachen/Documents/B_Cell_project/HVN037/cRel_RelA_correlation_sti/data_saved/CFP_nucleus_live_mask_value.mat')
load('/Users/yijiachen/Documents/B_Cell_project/HVN037/cRel_RelA_correlation_sti/data_saved/YFP_nucleus_live_mask_value.mat')

%% wc
%% 0 hour time point
% Initialize arrays to store the extracted values
% CFP_values = zeros(length(CFP_ConcIntensity_wc_tot), 1);
% YFP_values = zeros(length(YFP_ConcIntensity_wc_tot), 1);
%
% % Loop through each cell in the cell arrays
% for i = 1:length(CFP_ConcIntensity_wc_tot)
%     CFP_values(i) = CFP_ConcIntensity_wc_tot{i}(1,2);  % Extract value from first row, second column
%     YFP_values(i) = YFP_ConcIntensity_wc_tot{i}(1,2);  % Extract value from first row, second column
% end
CFP_values = double(CFP_nucleus_live_mask_all_t.nlm_medianint);
YFP_values = double(YFP_nucleus_live_mask_all_t.nlm_medianint);

correlation_plot = figure;
scatter(CFP_values, YFP_values, 50, 'filled','MarkerFaceColor','#009900');


p = polyfit(CFP_values, YFP_values, 1); % 1 for linear fit


% Generate x values for the line
x_fit = linspace(0, 100, 100);

% Evaluate the best-fit line at these x values
y_fit = polyval(p, x_fit);
% Calculate R^2 value
y_mean = mean(YFP_values);
SS_tot = sum((YFP_values - y_mean).^2);
SS_res = sum((YFP_values - polyval(p, CFP_values)).^2);
R_squared = 1 - (SS_res / SS_tot);

% Plot the best-fit line
hold on;
plot(x_fit, y_fit, 'r-', 'LineWidth', 2,'Color',[0.7,0.7,0.7]);
textString = sprintf('R^2: %.2f', R_squared);

% Update the legend
legend('Data', textString);
hold off;





% Add labels and title if needed
xlim([0 100])
ylim([0 150])

xlabel('CFP nuc MedianIntensity');
ylabel('YFP nuc MedianIntensity');
title('Correlation Plot of CFP vs YFP Nucleus Median Intensity with Stimulus');

saveas(correlation_plot,'/Users/yijiachen/Documents/B_Cell_project/HVN037/cRel_RelA_correlation_sti/images/CFP_YFP_correlation_plot/CFP_YFP_correlation_nuc_sti.png')


%
% %% 3 hour time point
% % Initialize arrays to store the extracted values
% CFP_values = zeros(length(CFP_ConcIntensity_wc_tot), 1);
% YFP_values = zeros(length(YFP_ConcIntensity_wc_tot), 1);
%
% % Loop through each cell in the cell arrays
% for i = 1:length(CFP_ConcIntensity_wc_tot)
%     CFP_values(i) = CFP_ConcIntensity_wc_tot{i}(19,2);  % Extract value from first row, second column
%     YFP_values(i) = YFP_ConcIntensity_wc_tot{i}(19,2);  % Extract value from first row, second column
% end
% correlation_3_hour = figure;
% scatter(CFP_values, YFP_values, 50, 'filled');
%
% % Add labels and title if needed
% xlim([2200 4800])
% ylim([500 1800])
%
% xlabel('CFP wc ConcIntensity');
% ylabel('YFP wc ConcIntensity');
% title('Correlation Plot of CFP vs YFP Wholecell Intensity at 3h');
%
% saveas(correlation_3_hour,'/Users/yijiachen/Documents/B_Cell_project/HVN037/images/CFP_YFP_correlation_plot/CFP_YFP_correlation_wc_3_hour.png')
%
% %% nuc
% %% 0 hour time point
% % Initialize arrays to store the extracted values
% CFP_values = zeros(length(CFP_MedianIntensity_nuc_tot), 1);
% YFP_values = zeros(length(YFP_MedianIntensity_nuc_tot), 1);
%
% % Loop through each cell in the cell arrays
% for i = 1:length(CFP_MedianIntensity_nuc_tot)
%     CFP_values(i) = CFP_MedianIntensity_nuc_tot{i}(1,2);  % Extract value from first row, second column
%     YFP_values(i) = YFP_MedianIntensity_nuc_tot{i}(1,2);  % Extract value from first row, second column
% end
%
% correlation_plot = figure;
% scatter(CFP_values, YFP_values, 50, 'filled');
%
% % Add labels and title if needed
% xlim([2200 4800])
% ylim([500 1800])
%
% xlabel('CFP nuc MedianIntensity');
% ylabel('YFP nuc MedianIntensity');
% title('Correlation Plot of CFP vs YFP Nucleus Intensity at 0h');
%
% saveas(correlation_plot,'/Users/yijiachen/Documents/B_Cell_project/HVN037/images/CFP_YFP_correlation_plot/CFP_YFP_correlation_nuc_0_hour.png')
%
%
%
% %% 3 hour time point
% % Initialize arrays to store the extracted values
% CFP_values = zeros(length(CFP_MedianIntensity_nuc_tot), 1);
% YFP_values = zeros(length(YFP_MedianIntensity_nuc_tot), 1);
%
% % Loop through each cell in the cell arrays
% for i = 1:length(CFP_MedianIntensity_nuc_tot)
%     CFP_values(i) = CFP_MedianIntensity_nuc_tot{i}(19,2);  % Extract value from first row, second column
%     YFP_values(i) = YFP_MedianIntensity_nuc_tot{i}(19,2);  % Extract value from first row, second column
% end
% correlation_3_hour = figure;
% scatter(CFP_values, YFP_values, 50, 'filled');
%
% % Add labels and title if needed
% xlim([2200 4800])
% ylim([500 1800])
%
% xlabel('CFP nuc MedianIntensity');
% ylabel('YFP nuc MedianIntensity');
% title('Correlation Plot of CFP vs YFP Nucleus Intensity at 3h');
%
% saveas(correlation_3_hour,'/Users/yijiachen/Documents/B_Cell_project/HVN037/images/CFP_YFP_correlation_plot/CFP_YFP_correlation_nuc_3_hour.png')
