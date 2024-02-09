load('/Users/yijiachen/Documents/B_Cell_project/HVN037/data_saved/HVN037_CFP_YFP_nuc_median_mean_int.mat')
load('/Users/yijiachen/Documents/B_Cell_project/HVN037/data_saved/HVN037_CFP_YFP_wc_median_mean_int.mat')

%% wc
%% 0 hour time point
% Initialize arrays to store the extracted values
CFP_values = zeros(length(CFP_ConcIntensity_wc_tot), 1);
YFP_values = zeros(length(YFP_ConcIntensity_wc_tot), 1);

% Loop through each cell in the cell arrays
for i = 1:length(CFP_ConcIntensity_wc_tot)
    CFP_values(i) = CFP_ConcIntensity_wc_tot{i}(1,2);  % Extract value from first row, second column
    YFP_values(i) = YFP_ConcIntensity_wc_tot{i}(1,2);  % Extract value from first row, second column
end

correlation_0_hour = figure;
scatter(CFP_values, YFP_values, 50, 'filled','MarkerFaceColor','#009900'); 




% Perform the linear fitting with the cleaned data
p = polyfit(CFP_values, YFP_values, 1); % 1 for linear fit

% Generate x values for the line
x_fit = linspace(2500, 5000, 100);

% Evaluate the best-fit line at these x values
y_fit = polyval(p, x_fit);

% Calculate R^2 value using the cleaned data
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
xlim([2500 5000])
ylim([500 2000])

xlabel('CFP wc ConcIntensity');
ylabel('YFP wc ConcIntensity');
title('Correlation Plot of CFP vs YFP Wholecell Intensity at 0h');

saveas(correlation_0_hour,'/Users/yijiachen/Documents/B_Cell_project/HVN037/images/CFP_YFP_correlation_plot/CFP_YFP_correlation_wc_0_hour.png')



%% 3 hour time point
% Initialize arrays to store the extracted values
CFP_values = zeros(length(CFP_ConcIntensity_wc_tot), 1);
YFP_values = zeros(length(YFP_ConcIntensity_wc_tot), 1);

% Loop through each cell in the cell arrays
for i = 1:length(CFP_ConcIntensity_wc_tot)
    CFP_values(i) = CFP_ConcIntensity_wc_tot{i}(19,2);  % Extract value from first row, second column
    YFP_values(i) = YFP_ConcIntensity_wc_tot{i}(19,2);  % Extract value from first row, second column
end
correlation_3_hour = figure;
scatter(CFP_values, YFP_values, 50, 'filled','MarkerFaceColor','#009900'); 

p = polyfit(CFP_values, YFP_values, 1); % 1 for linear fit


% Generate x values for the line
x_fit = linspace(min(CFP_values), 5000, 100);

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
xlim([2500 5000])
ylim([500 2000])

xlabel('CFP wc ConcIntensity');
ylabel('YFP wc ConcIntensity');
title('Correlation Plot of CFP vs YFP Wholecell Intensity at 3h');

saveas(correlation_3_hour,'/Users/yijiachen/Documents/B_Cell_project/HVN037/images/CFP_YFP_correlation_plot/CFP_YFP_correlation_wc_3_hour.png')


%% nuc
%% 0 hour time point
% Initialize arrays to store the extracted values
CFP_values = zeros(length(CFP_MedianIntensity_nuc_tot), 1);
YFP_values = zeros(length(YFP_MedianIntensity_nuc_tot), 1);

% Loop through each cell in the cell arrays
for i = 1:length(CFP_MedianIntensity_nuc_tot)
    CFP_values(i) = CFP_MedianIntensity_nuc_tot{i}(1,2);  % Extract value from first row, second column
    YFP_values(i) = YFP_MedianIntensity_nuc_tot{i}(1,2);  % Extract value from first row, second column
end

correlation_0_hour = figure;
scatter(CFP_values, YFP_values, 50, 'filled','MarkerFaceColor','#009900'); 

validIndices = ~isnan(CFP_values) & ~isnan(YFP_values);
cleaned_CFP_values = CFP_values(validIndices);
cleaned_YFP_values = YFP_values(validIndices);


% Perform the linear fitting with the cleaned data
p = polyfit(cleaned_CFP_values, cleaned_YFP_values, 1); % 1 for linear fit

% Generate x values for the line
x_fit = linspace(2500, 5000, 100);

% Evaluate the best-fit line at these x values
y_fit = polyval(p, x_fit);

% Calculate R^2 value using the cleaned data
y_mean = mean(cleaned_YFP_values);
SS_tot = sum((cleaned_YFP_values - y_mean).^2);
SS_res = sum((cleaned_YFP_values - polyval(p, cleaned_CFP_values)).^2);
R_squared = 1 - (SS_res / SS_tot);

% Plot the best-fit line
hold on;
plot(x_fit, y_fit, 'r-', 'LineWidth', 2,'Color',[0.7,0.7,0.7]);
textString = sprintf('R^2: %.2f', R_squared);

% Update the legend
legend('Data', textString);
hold off;

% Add labels and title if needed
xlim([2500 5000])
ylim([500 2000])

xlabel('CFP nuc MedianIntensity');
ylabel('YFP nuc MedianIntensity');
title('Correlation Plot of CFP vs YFP Nucleus Intensity at 0h');

saveas(correlation_0_hour,'/Users/yijiachen/Documents/B_Cell_project/HVN037/images/CFP_YFP_correlation_plot/CFP_YFP_correlation_nuc_0_hour.png')



%% 3 hour time point
% Initialize arrays to store the extracted values
CFP_values = zeros(length(CFP_MedianIntensity_nuc_tot), 1);
YFP_values = zeros(length(YFP_MedianIntensity_nuc_tot), 1);

% Loop through each cell in the cell arrays
for i = 1:length(CFP_MedianIntensity_nuc_tot)
    CFP_values(i) = CFP_MedianIntensity_nuc_tot{i}(19,2);  % Extract value from first row, second column
    YFP_values(i) = YFP_MedianIntensity_nuc_tot{i}(19,2);  % Extract value from first row, second column
end
correlation_3_hour = figure;
scatter(CFP_values, YFP_values, 50, 'filled','MarkerFaceColor','#009900'); 

validIndices = ~isnan(CFP_values) & ~isnan(YFP_values);
cleaned_CFP_values = CFP_values(validIndices);
cleaned_YFP_values = YFP_values(validIndices);

% Perform the linear fitting with the cleaned data
p = polyfit(cleaned_CFP_values, cleaned_YFP_values, 1); % 1 for linear fit

% Generate x values for the line
x_fit = linspace(2500, 5000, 100);

% Evaluate the best-fit line at these x values
y_fit = polyval(p, x_fit);

% Calculate R^2 value using the cleaned data
y_mean = mean(cleaned_YFP_values);
SS_tot = sum((cleaned_YFP_values - y_mean).^2);
SS_res = sum((cleaned_YFP_values - polyval(p, cleaned_CFP_values)).^2);
R_squared = 1 - (SS_res / SS_tot);

% Plot the best-fit line
hold on;
plot(x_fit, y_fit, 'r-', 'LineWidth', 2,'Color',[0.7,0.7,0.7]);
textString = sprintf('R^2: %.2f', R_squared);

% Update the legend
legend('Data', textString);
hold off;

% Add labels and title if needed
xlim([2500 5000])
ylim([500 2000])

xlabel('CFP nuc MedianIntensity');
ylabel('YFP nuc MedianIntensity');
title('Correlation Plot of CFP vs YFP Nucleus Intensity at 3h');

saveas(correlation_3_hour,'/Users/yijiachen/Documents/B_Cell_project/HVN037/images/CFP_YFP_correlation_plot/CFP_YFP_correlation_nuc_3_hour.png')
