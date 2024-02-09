load('/Users/yijiachen/Documents/B_Cell_project/HVN037/data_saved/HVN037_CFP_YFP_nuc_median_mean_int.mat')
load('/Users/yijiachen/Documents/B_Cell_project/HVN037/data_saved/HVN037_CFP_YFP_wc_median_mean_int.mat')

%% CFP nuc median int traj
CFP_nuclear_medianint_plot = figure;
for s = 1
    for cell = [5,7,9,17,24,25,26]
        %% CFP nuclear median intensity traj

        valid_indices_CFP = ~isnan(CFP_MedianIntensity_nuc_tot{cell,1}(:, 2)) & CFP_MedianIntensity_nuc_tot{cell,1}(:, 2) ~= 0;
        plot(CFP_MedianIntensity_nuc_tot{cell,1}(valid_indices_CFP, 1), CFP_MedianIntensity_nuc_tot{cell,1}(valid_indices_CFP, 2), 'LineWidth', 2);
        hold on
    end
    
    CFP_rgb_colors = [0,204,255; 0,0,255; 0,102,255; 0,153,255; 153,204,255; 51,51,204; 51,102,204];
    CFP_colors = CFP_rgb_colors / 255;
    ax = gca;
    ax.ColorOrder = CFP_colors;

    ylim([2500 5000])
    xlim([0 360])
    xticks([1 60:60:360])
    xticklabels({'0', '1', '2', '3', '4', '5', '6'});
    title('CFP nuclear median intensity')
    xlabel('time (h)')
    ylabel('CFP nuclear median intensity')
end
saveas(CFP_nuclear_medianint_plot,'/Users/yijiachen/Documents/B_Cell_project/HVN037/images/CFP_nuc_traj/fluor_plot_medianint.png')

%% CFP wc concentration int traj
CFP_wc_concint_plot = figure;
for s = 1
    for cell = [5,7,9,17,24,25,26]
        %% CFP wc concentration intensity traj

        valid_indices_CFP = ~isnan(CFP_ConcIntensity_wc_tot{cell,1}(:, 2)) & CFP_ConcIntensity_wc_tot{cell,1}(:, 2) ~= 0;

        plot(CFP_ConcIntensity_wc_tot{cell,1}(valid_indices_CFP, 1), CFP_ConcIntensity_wc_tot{cell,1}(valid_indices_CFP, 2), 'LineWidth', 2);

        hold on

    end

    CFP_rgb_colors = [0,204,255; 0,0,255; 0,102,255; 0,153,255; 153,204,255; 51,51,204; 51,102,204];
    CFP_colors = CFP_rgb_colors / 255;
    ax = gca;
    ax.ColorOrder = CFP_colors;

    ylim([2500 5000])
    xlim([0 360])
    xticks([1 60:60:360])
    % %     %     xticklabels({'0','2','4','6','8','10','12','14','16'})
    % xticklabels({'0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16'});
    xticklabels({'0', '1', '2', '3', '4', '5', '6'});

    title('CFP wc concentration intensity')
    xlabel('time (h)')
    ylabel('CFP wc concentration intensity')
    % legend ('CFP','YFP')
end
saveas(CFP_wc_concint_plot,'/Users/yijiachen/Documents/B_Cell_project/HVN037/images/CFP_wc_traj/fluor_plot_concint.png')


%% YFP nuc median int traj
YFP_nuclear_medianint_plot = figure;
for s = 1
    for cell = [5,7,9,17,24,25,26]
        %% YFP nuclear median intensity traj

        valid_indices_YFP = ~isnan(YFP_MedianIntensity_nuc_tot{cell,1}(:, 2)) & YFP_MedianIntensity_nuc_tot{cell,1}(:, 2) ~= 0;

        plot(YFP_MedianIntensity_nuc_tot{cell,1}(valid_indices_YFP, 1), YFP_MedianIntensity_nuc_tot{cell,1}(valid_indices_YFP, 2), 'LineWidth', 2);

        hold on

    end

    YFP_rgb_colors = [255,204,0; 255,255,0; 204,153,0; 255,153,0; 255,255,102; 184,180,0; 125,122,0];
    YFP_colors = YFP_rgb_colors / 255;
    ax = gca;
    ax.ColorOrder = YFP_colors;

    ylim([500 2000])
    xlim([0 360])
    xticks([1 60:60:360])
    % %     %     xticklabels({'0','2','4','6','8','10','12','14','16'})
    % xticklabels({'0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16'});
    xticklabels({'0', '1', '2', '3', '4', '5', '6'});

    title('YFP nuclear median intensity')
    xlabel('time (h)')
    ylabel('YFP nuclear median intensity')
    % legend ('YFP','YFP')
end
saveas(YFP_nuclear_medianint_plot,'/Users/yijiachen/Documents/B_Cell_project/HVN037/images/YFP_nuc_traj/fluor_plot_medianint.png')

%% YFP wc concentration int traj
YFP_wc_concint_plot = figure;
for s = 1
    for cell = [5,7,9,17,24,25,26]
        %% YFP wc concentration intensity traj

        valid_indices_YFP = ~isnan(YFP_ConcIntensity_wc_tot{cell,1}(:, 2)) & YFP_ConcIntensity_wc_tot{cell,1}(:, 2) ~= 0;

        plot(YFP_ConcIntensity_wc_tot{cell,1}(valid_indices_YFP, 1), YFP_ConcIntensity_wc_tot{cell,1}(valid_indices_YFP, 2), 'LineWidth', 2);

        hold on

    end

    YFP_rgb_colors = [255,204,0; 255,255,0; 204,153,0; 255,153,0; 255,255,102; 184,180,0; 125,122,0];
    YFP_colors = YFP_rgb_colors / 255;

    ax = gca;
    ax.ColorOrder = YFP_colors;

    ylim([500 2000])
    xlim([0 360])
    xticks([1 60:60:360])
    % %     %     xticklabels({'0','2','4','6','8','10','12','14','16'})
    % xticklabels({'0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16'});
    xticklabels({'0', '1', '2', '3', '4', '5', '6'});

    title('YFP wc concentration intensity')
    xlabel('time (h)')
    ylabel('YFP wc concentration intensity')
    % legend ('YFP','YFP')
end
saveas(YFP_wc_concint_plot,'/Users/yijiachen/Documents/B_Cell_project/HVN037/images/YFP_wc_traj/fluor_plot_concint.png')

