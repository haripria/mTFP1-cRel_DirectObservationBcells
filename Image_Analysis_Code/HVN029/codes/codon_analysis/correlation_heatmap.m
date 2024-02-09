% Load Data
load ('/Users/yijiachen/Documents/B_Cell_project/HVN029/data_saved/HVN029_CFP_nuc_peak_amp_time.mat')
load ('/Users/yijiachen/Documents/B_Cell_project/HVN029/data_saved/HVN029_YFP_nuc_peak_amp_time.mat')

% Define Variables
CFP_YFP_variables = [CFP_peak_amp, CFP_peak_time, CFP_basal_amp, CFP_fold_change, CFP_pb_diff, YFP_peak_amp, YFP_peak_time, YFP_basal_amp];
% CFP_YFP_variables = CFP_YFP_variables(CFP_YFP_variables(:, 4) < 4, :);
% Variable Names
CFP_YFP_variable_names = {'cRel Peak Amp', 'cRel Peak Time', 'cRel Basal Amp', 'cRel Fold Change', 'cRel Peak-Basal Diff', 'cMyc Peak Amp', 'cMyc Peak Time', 'cMyc Basal Amp'};

% Calculate Correlations and p-values
[CFP_YFP_r, CFP_YFP_pval] = corr(CFP_YFP_variables, 'Rows','complete');

% Create a figure for the heatmap
CFP_YFP_corr_heatmap = figure;

% Create the heatmap manually
imagesc(CFP_YFP_r);
colorbar;
caxis([-1,1])
colormap("turbo")
axis square;
set(gca, 'XTick', 1:length(CFP_YFP_variable_names), 'XTickLabel', CFP_YFP_variable_names, 'YTick', 1:length(CFP_YFP_variable_names), 'YTickLabel', CFP_YFP_variable_names);
xtickangle(45);
title('CFP YFP Correlation Coefficients WT');

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
saveas(CFP_YFP_corr_heatmap, '/Users/yijiachen/Documents/B_Cell_project/HVN029/images/CFP_YFP_codon_correlation_heatmap.png');

% saveas(CFP_YFP_corr_heatmap, '/Users/yijiachen/Documents/B_Cell_project/HVN029/images/CFP_YFP_codon_correlation_heatmap_filtered.png');
