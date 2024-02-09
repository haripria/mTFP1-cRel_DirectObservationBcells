% Load Data
load ('/Users/yijiachen/Documents/B_Cell_project/HVN029/data_saved/HVN029_CFP_nuc_peak_amp_time.mat')
load ('/Users/yijiachen/Documents/B_Cell_project/HVN029/data_saved/HVN029_YFP_nuc_peak_amp_time.mat')
CFP_basal_amp_wt = CFP_basal_amp;
CFP_fold_change_wt = CFP_fold_change;
CFP_peak_amp_wt = CFP_peak_amp;
CFP_peak_time_wt = CFP_peak_time;
clear CFP_basal_amp CFP_fold_change CFP_peak_amp CFP_peak_time

YFP_basal_amp_wt = YFP_basal_amp;
YFP_peak_amp_wt = YFP_peak_amp;
YFP_peak_time_wt = YFP_peak_time;
clear YFP_basal_amp YFP_peak_amp YFP_peak_time

load ('/Users/yijiachen/Documents/B_Cell_project/HVN032/data_saved/HVN032_CFP_nuc_peak_amp_time.mat')
load ('/Users/yijiachen/Documents/B_Cell_project/HVN032/data_saved/HVN032_YFP_nuc_peak_amp_time.mat')
CFP_basal_amp_eko = CFP_basal_amp;
CFP_fold_change_eko = CFP_fold_change;
CFP_peak_amp_eko = CFP_peak_amp;
CFP_peak_time_eko = CFP_peak_time;
clear CFP_basal_amp CFP_fold_change CFP_peak_amp CFP_peak_time

YFP_basal_amp_eko = YFP_basal_amp;
YFP_peak_amp_eko = YFP_peak_amp;
YFP_peak_time_eko = YFP_peak_time;
clear YFP_basal_amp YFP_peak_amp YFP_peak_time

% Define Variables
CFP_basal_amp = [CFP_basal_amp_wt; CFP_basal_amp_eko];
CFP_fold_change = [CFP_fold_change_wt; CFP_fold_change_eko];
CFP_peak_amp = [CFP_peak_amp_wt; CFP_peak_amp_eko];
CFP_peak_time = [CFP_peak_time_wt; CFP_peak_time_eko];
YFP_basal_amp = [YFP_basal_amp_wt; YFP_basal_amp_eko];
YFP_peak_amp = [YFP_peak_amp_wt; YFP_peak_amp_eko];
YFP_peak_time = [YFP_peak_time_wt; YFP_peak_time_eko];


% Define Variables
CFP_YFP_variables = [CFP_peak_amp, CFP_peak_time, CFP_basal_amp, CFP_fold_change, YFP_peak_amp, YFP_peak_time, YFP_basal_amp];
CFP_YFP_variables = CFP_YFP_variables(CFP_YFP_variables(:, 4) < 4, :);
% Variable Names
CFP_YFP_variable_names = {'cRel Peak Amp', 'cRel Peak Time', 'cRel Basal Amp', 'cRel Fold Change', 'cMyc Peak Amp', 'cMyc Peak Time', 'cMyc Basal Amp'};

% Calculate Correlations and p-values
[CFP_YFP_r, CFP_YFP_pval] = corr(CFP_YFP_variables, 'Rows','complete');

% Create a figure for the heatmap
CFP_YFP_corr_heatmap = figure;

% Create the heatmap manually
imagesc(CFP_YFP_r);
colorbar;
axis square;
set(gca, 'XTick', 1:length(CFP_YFP_variable_names), 'XTickLabel', CFP_YFP_variable_names, 'YTick', 1:length(CFP_YFP_variable_names), 'YTickLabel', CFP_YFP_variable_names);
xtickangle(45);
title('CFP YFP Correlation Coefficients WT eKO together');

% Overlay r values and stars for significant p-values
[rows, cols] = size(CFP_YFP_r);
for i = 1:rows
    for j = 1:cols
%         textValue = sprintf('%.2f', CFP_YFP_r(i, j));
        rValue = sprintf('%.2f', CFP_YFP_r(i, j));
        pValue = sprintf('%.3f', CFP_YFP_pval(i, j));
        textValue = [rValue, newline, pValue]; % Combine r and p values in one string

        if CFP_YFP_pval(i, j) < 0.05
            textValue = [textValue, '*']; % Append star if p-value is significant
        end
        text(j, i, textValue, 'HorizontalAlignment', 'Center', 'Color', 'white');
    end
end

% Save the heatmap
saveas(CFP_YFP_corr_heatmap,'/Users/yijiachen/Documents/B_Cell_project/images/CFP_YFP_codon_correlation_heatmap_wt_eko_filtered.png')


% 
% 
% 
% % Calculate Correlations and p-values for CFP and YFP
% [CFP_YFP_r, CFP_YFP_pval] = corr(CFP_YFP_variables, 'Rows','complete');
% % [YFP_r, YFP_pval] = corr(YFP_variables, 'Rows','complete');
% 
% % Plot Heatmap for Correlation Coefficients and p-values for CFP
% CFP_YFP_corr_heatmap = figure;
% subplot(1,2,1);
% heatmap(CFP_YFP_variable_names, CFP_YFP_variable_names, CFP_YFP_r);
% title('CFP YFP Correlation Coefficients WT eKO together');
% 
% subplot(1,2,2);
% heatmap(CFP_YFP_variable_names, CFP_YFP_variable_names, CFP_YFP_pval);
% title('CFP YFP p-values');
% 
% saveas(CFP_YFP_corr_heatmap,'/Users/yijiachen/Documents/B_Cell_project/images/CFP_YFP_codon_correlation_heatmap_wt_eko.png')
