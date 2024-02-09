load('/Users/yijiachen/Documents/B_Cell_project/HVN032/data_saved/HVN032_CFP_YFP_nuc_median_mean_int.mat')
load('/Users/yijiachen/Documents/B_Cell_project/HVN032/data_saved/HVN032_s1_s4_to_s9_tot_cell.mat')

tot_cell_num = 1;
CFP_bg = 2200;

filtered_cell = [29,37,45,48];
for filtered_cell_num = 1:size(filtered_cell,2)
    tot_cell(tot_cell == filtered_cell(filtered_cell_num)) = [];
end

%% visualization
CFP_nuclear_medianint_plot_smoothed = figure;

%% analysis
for cell = tot_cell
    %[1,2,3,4,6,7,8,9,15]
    %% CFP nuclear median intensity traj


  
    valid_indices_CFP = ~isnan(CFP_MedianIntensity_nuc_tot{cell,1}(:, 2)) & CFP_MedianIntensity_nuc_tot{cell,1}(:, 2) ~= 0 & CFP_MedianIntensity_nuc_tot{cell,1}(:, 1) <= 961;

    valid_time_CFP = CFP_MedianIntensity_nuc_tot{cell,1}(valid_indices_CFP, 1);
    valid_medianint_CFP = CFP_MedianIntensity_nuc_tot{cell,1}(valid_indices_CFP, 2);

    %% bg rescale
    % % %     valid_medianint_CFP = CFP_MedianIntensity_nuc_tot{cell,1}(valid_indices_CFP, 2);
    % %     valid_medianint_CFP_bg_rescale = (valid_medianint_CFP-CFP_bg)./1000; %% set bg to 0, every 1000 increase scale to 1
    % % %     plot(CFP_MedianIntensity_nuc_tot{cell,1}(valid_indices_CFP, 1), valid_medianint_CFP_bg_rescale, 'LineWidth', 2);

    valid_medianint_CFP_bg_rescale = (valid_medianint_CFP-CFP_bg); %% set bg to 0


    valid_time_int_CFP = [valid_time_CFP valid_medianint_CFP_bg_rescale];
    order = 3;
    framelen = 11;

    CFP_sgf = sgolayfilt(valid_time_int_CFP,order,framelen);

    CFP_sgf_tot{tot_cell_num} = CFP_sgf;

    % %      [minValue, minIndex] = min(valid_time_CFP);
    % %     CFP_basal_amp(tot_cell_num,1) = valid_medianint_CFP_bg_rescale(minIndex, 1);

    [minValue, minIndex] = min(valid_time_CFP(valid_time_CFP>=21));
    CFP_basal_time(tot_cell_num,1)  = valid_time_CFP(valid_time_CFP==minValue);
    CFP_basal_amp(tot_cell_num,1)  =CFP_sgf(valid_time_CFP==minValue, 2);

    %% visualization
    plot(valid_time_int_CFP(:,1),valid_time_int_CFP(:,2),':','LineWidth',2)
    hold on
    plot(CFP_sgf(:,1),CFP_sgf(:,2),'.-','LineWidth',2)
    hold on

    tot_cell_num = tot_cell_num +1;
end

CFP_peak_amp = zeros(tot_cell_num-1,1);
CFP_peak_time = zeros(tot_cell_num-1,1);

%% extract the peak(s) amp and time
for c_num = 1: (tot_cell_num-1)
    CFP_sgf_time = CFP_sgf_tot{c_num}(:,1);
    CFP_sgf_int = CFP_sgf_tot{c_num}(:,2);
    [CFP_pks,locs] = findpeaks(CFP_sgf_int(CFP_sgf_time<602));
    [CFP_max_peak, CFP_max_idx] = max(CFP_pks); % Find the highest peak
    CFP_peak_location = locs(CFP_max_idx);
    CFP_peak_amp(c_num,1) = CFP_sgf_tot{c_num}(CFP_peak_location,2);
    CFP_peak_time(c_num,1) = CFP_sgf_tot{c_num}(CFP_peak_location,1);
end


CFP_fold_change = CFP_peak_amp./CFP_basal_amp;
CFP_pb_diff = CFP_peak_amp - CFP_basal_amp;

save('/Users/yijiachen/Documents/B_Cell_project/HVN032/data_saved/HVN032_CFP_nuc_peak_amp_time.mat','CFP_peak_amp','CFP_peak_time','CFP_basal_amp','CFP_basal_time','CFP_fold_change','CFP_pb_diff')


%% visualization
scatter(CFP_peak_time,CFP_peak_amp,100,'filled')
scatter(CFP_basal_time,CFP_basal_amp,100,'filled')

CFP_rgb_colors = [0,204,255; 0,0,255; 0,102,255; 0,153,255; 153,204,255; 51,51,204; 51,102,204; 0,255,255; 102,255,255];
CFP_colors = CFP_rgb_colors / 255;
box off; % Turn off the box surrounding the plot
ax = gca;
ax.ColorOrder = CFP_colors;
ax.TickDir = 'in'; % Set the direction of the ticks

% ylim([2000 3200])
% ylim([0 1])
xlim([0 960])
xticks([1 60:60:960])
xticklabels({'0', '1', '2', '3', '4', '5', '6','7', '8', '9', '10', '11', '12', '13','14', '15', '16'});
title('CFP nuclear median intensity')
xlabel('time (h)')
ylabel('CFP nuclear median intensity')
legend('signal','sgolay')
saveas(CFP_nuclear_medianint_plot_smoothed,'/Users/yijiachen/Documents/B_Cell_project/HVN032/images/CFP_nuc_rep_traj/fluor_plot_medianint_smoothed_scalednew.png')
% % % % saveas(CFP_nuclear_medianint_plot_smoothed,sprintf('/Users/yijiachen/Documents/B_Cell_project/HVN032/images/CFP_nuc_rep_traj/fluor_plot_medianint_smoothed_c%d.png',cell))




