% Load Data
load ('/Users/yijiachen/Documents/B_Cell_project/HVN029/data_saved/HVN029_CFP_nuc_peak_amp_time.mat')
load ('/Users/yijiachen/Documents/B_Cell_project/HVN029/data_saved/HVN029_YFP_nuc_peak_amp_time.mat')

% Define Variables
CFP_YFP_variables = [CFP_peak_amp, CFP_peak_time, CFP_basal_amp, CFP_fold_change, CFP_pb_diff, YFP_peak_amp, YFP_peak_time, YFP_basal_amp];

CFP_YFP_variables = CFP_YFP_variables(CFP_YFP_variables(:, 4) < 4, :);
% Variable Names
CFP_YFP_variable_names = {'cRel Peak Amp', 'cRel Peak Time', 'cRel Basal Amp', 'cRel Fold Change', 'cRel Peak-Basal Diff','cMyc Peak Amp', 'cMyc Peak Time', 'cMyc Basal Amp'};

% Calculate Correlations and p-values
[CFP_YFP_r, CFP_YFP_pval] = corr(CFP_YFP_variables, 'Rows','complete');

% Number of Variables
numVars = size(CFP_YFP_variables, 2);

% Loop through pairs of variables
for i = 1:numVars
    for j = 1:i-1 % Only consider lower half to avoid redundant pairs
        %         if CFP_YFP_pval(i, j) < 0.05
        % Create a new figure for each significant pair
        scatter_plot = figure;

        % Scatter Plot
        h = scatter(CFP_YFP_variables(:,j), CFP_YFP_variables(:,i),50, 'black','filled');



        % Set the alpha transparency to 0.5
        h.MarkerFaceAlpha = 0.7;
        h.MarkerEdgeAlpha = 0.7;

        xlabel(CFP_YFP_variable_names{j});
        ylabel(CFP_YFP_variable_names{i});
        title(sprintf('Scatter Plot of %s vs %s', CFP_YFP_variable_names{j}, CFP_YFP_variable_names{i}));

        % Best Fit Line
        coeffs = polyfit(CFP_YFP_variables(:,j), CFP_YFP_variables(:,i), 1); % Linear fit

        xlims = xlim;

        %         fittedX = linspace(min(CFP_YFP_variables(:,j)), max(CFP_YFP_variables(:,j)), 200);
        fittedX = linspace(xlims(1), xlims(2), 200);
        fittedY = polyval(coeffs, fittedX);
        hold on;
        plot(fittedX, fittedY, 'r','LineWidth',2);
        hold off;

        % Display Correlation Coefficient
        annotation('textbox', [.15 .8 .1 .1], 'String', sprintf('r = %.2f, p = %.2f', CFP_YFP_r(i, j), CFP_YFP_pval(i,j)), 'EdgeColor', 'none');
        saveas(scatter_plot, sprintf('/Users/yijiachen/Documents/B_Cell_project/HVN029/images/correlation_scatter_filtered/%s_%s.png', CFP_YFP_variable_names{j}, CFP_YFP_variable_names{i}));
        %         end
    end
end
