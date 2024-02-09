load('/Users/yijiachen/Documents/B_Cell_project/HVN029/data_saved/HVN029_CFP_YFP_nuc_median_mean_int.mat')
load('/Users/yijiachen/Documents/B_Cell_project/HVN029/data_saved/HVN029_s2_to_s8_tot_cell.mat')
load('/Users/yijiachen/Documents/B_Cell_project/HVN029/data_saved/HVN029_CFP_nuc_peak_amp_time.mat')
load('/Users/yijiachen/Documents/B_Cell_project/HVN029/data_saved/HVN029_YFP_nuc_peak_amp_time.mat')


tot_cell_num = 1;
CFP_bg = 2300;
YFP_bg = 700;
filtered_cell = [85,72,62,54,49,42,38,10,8,7,6];
for filtered_cell_num = 1:size(filtered_cell,2)
    tot_cell(tot_cell == filtered_cell(filtered_cell_num)) = [];
end


%% CFP nuc median int traj

for cell = tot_cell
    CFP_YFP_nuclear_medianint_plot = figure;
    %% CFP nuclear median intensity traj
    %     valid_indices_CFP = ~isnan(CFP_MedianIntensity_nuc_tot{cell,1}(:, 2)) & CFP_MedianIntensity_nuc_tot{cell,1}(:, 2) ~= 0;
    valid_indices_CFP = ~isnan(CFP_MedianIntensity_nuc_tot{cell,1}(:, 2)) & CFP_MedianIntensity_nuc_tot{cell,1}(:, 2) ~= 0 & CFP_MedianIntensity_nuc_tot{cell,1}(:, 1) <= 961;

    valid_time_CFP = CFP_MedianIntensity_nuc_tot{cell,1}(valid_indices_CFP, 1);
    %% bg rescale
    valid_medianint_CFP = CFP_MedianIntensity_nuc_tot{cell,1}(valid_indices_CFP, 2);
    %     valid_medianint_CFP_bg_rescale = (valid_medianint_CFP-CFP_bg)./1000; %% set bg to 0, every 1000 increase scale to 1
    %     plot(CFP_MedianIntensity_nuc_tot{cell,1}(valid_indices_CFP, 1), valid_medianint_CFP_bg_rescale, 'LineWidth', 2);
    %     hold on

    valid_medianint_CFP_bg_rescale = (valid_medianint_CFP-CFP_bg); %% set bg to 0

    valid_time_int_CFP = [valid_time_CFP valid_medianint_CFP_bg_rescale];
    order = 3;
    framelen = 11;

    CFP_sgf = sgolayfilt(valid_time_int_CFP,order,framelen);
    CFP_sgf_tot{tot_cell_num} = CFP_sgf;


    %% visualization
    % Create the first plot with the first y-axis
    yyaxis left
    plot(valid_time_int_CFP(:,1),valid_time_int_CFP(:,2),':','LineWidth',2)
    hold on
    plot(CFP_sgf(:,1),CFP_sgf(:,2),'.-','LineWidth',2)
    hold on
    scatter(CFP_peak_time(tot_cell_num),CFP_peak_amp(tot_cell_num),100,'filled')
    scatter(CFP_basal_time(tot_cell_num),CFP_basal_amp(tot_cell_num),100,'filled')

    %     ylim(yyaxis left, [ymin1 ymax1])
    %     ylim([0 3200-CFP_bg])
    ylabel('Y-axis for CFP')

    CFP_rgb_colors = [0,204,255];
    CFP_colors = CFP_rgb_colors / 255;
    box off; % Turn off the box surrounding the plot
    ax = gca;
    ax.ColorOrder = CFP_colors;
    ax.TickDir = 'in'; % Set the direction of the ticks
        ax.YAxis(1).Color = CFP_colors; 

    %% YFP nuclear median intensity traj
    valid_indices_YFP = ~isnan(YFP_MedianIntensity_nuc_tot{cell,1}(:, 2)) & YFP_MedianIntensity_nuc_tot{cell,1}(:, 2) ~= 0 & YFP_MedianIntensity_nuc_tot{cell,1}(:, 1) <= 961;
    valid_time_YFP = YFP_MedianIntensity_nuc_tot{cell,1}(valid_indices_YFP, 1);
    %     valid_medianint_YFP = YFP_MedianIntensity_nuc_tot{cell,1}(valid_indices_YFP, 2);
    %% bg rescale
    valid_medianint_YFP = YFP_MedianIntensity_nuc_tot{cell,1}(valid_indices_YFP, 2);
    valid_medianint_YFP_bg_rescale = (valid_medianint_YFP-YFP_bg); %% set bg to 0
    valid_time_int_YFP = [valid_time_YFP valid_medianint_YFP_bg_rescale];
    order = 3;
    framelen = 11;
    YFP_sgf = sgolayfilt(valid_time_int_YFP,order,framelen);
    YFP_sgf_tot{tot_cell_num} = YFP_sgf;


    %% visualization
    yyaxis right
    plot(valid_time_int_YFP(:,1),valid_time_int_YFP(:,2),':','LineWidth',2)
    hold on
    plot(YFP_sgf(:,1),YFP_sgf(:,2),'.-','LineWidth',2)
    hold on

    scatter(YFP_peak_time(tot_cell_num),YFP_peak_amp(tot_cell_num),100,'filled')
    scatter(YFP_basal_time(tot_cell_num),YFP_basal_amp(tot_cell_num),100,'filled')
    ylabel('Y-axis for YFP')

 YFP_rgb_colors = [255,204,0];
    YFP_colors = YFP_rgb_colors / 255;
    box off; % Turn off the box surrounding the plot
    ax = gca;
    ax.ColorOrder = YFP_colors;
    ax.TickDir = 'in'; % Set the direction of the ticks
    ax.YAxis(2).Color = YFP_colors; 





    %     plot(YFP_MedianIntensity_nuc_tot{cell,1}(valid_indices_YFP, 1), YFP_MedianIntensity_nuc_tot{cell,1}(valid_indices_YFP, 2), 'LineWidth', 2);

    %     hold on

%     CFP_rgb_colors = [0,204,255; 0,204,255; 0,204,255; 0,204,255; 255,204,0; 255,204,0;255,204,0; 255,204,0];
%     CFP_colors = CFP_rgb_colors / 255;
%     box off; % Turn off the box surrounding the plot
%     ax = gca;
%     ax.ColorOrder = CFP_colors;
%     ax.TickDir = 'in'; % Set the direction of the ticks

    %     % ylim([2000 3200])
    %     ylim([0 1])





    % Optionally set axis limits and other properties
    % xlim([xmin xmax])
    %
    % ylim(yyaxis right, [ymin2 ymax2])


    xlim([0 960])
    xticks([1 60:60:960])
    xticklabels({'0', '1', '2', '3', '4', '5', '6','7', '8', '9', '10', '11', '12', '13','14', '15', '16'});
    title(sprintf('CFP/YFP nuclear median intensity c%d',tot_cell(tot_cell_num)))
    xlabel('time (h)')
%     ylabel('CFP/YFP nuclear median intensity')



    saveas(CFP_YFP_nuclear_medianint_plot,sprintf('/Users/yijiachen/Documents/B_Cell_project/HVN029/images/CFP_YFP_nuc_rep_traj_two_yaxis/fluor_plot_medianint_c%d.png',tot_cell(tot_cell_num)))
    tot_cell_num=tot_cell_num+1;
end

