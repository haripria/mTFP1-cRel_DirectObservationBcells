load('/Users/yijiachen/Documents/B_Cell_project/HVN029/data_saved/HVN029_CFP_YFP_nuc_median_mean_int.mat')
load('/Users/yijiachen/Documents/B_Cell_project/HVN029/data_saved/HVN029_s2_to_s8_tot_cell.mat')

tot_cell_num = 1;
YFP_bg = 700;

filtered_cell = [85,72,62,54,49,42,38,10,8,7,6];
for filtered_cell_num = 1:size(filtered_cell,2)
    tot_cell(tot_cell == filtered_cell(filtered_cell_num)) = [];
end

%% visualization
YFP_nuclear_medianint_plot_smoothed = figure;

%% analysis
for cell = tot_cell
    %[6,7,9,10,23,27,28,31,32,33]
    %% YFP nuclear median intensity traj

    valid_indices_YFP = ~isnan(YFP_MedianIntensity_nuc_tot{cell,1}(:, 2)) & YFP_MedianIntensity_nuc_tot{cell,1}(:, 2) ~= 0 & YFP_MedianIntensity_nuc_tot{cell,1}(:, 1) <= 961;

    valid_time_YFP = YFP_MedianIntensity_nuc_tot{cell,1}(valid_indices_YFP, 1);
    valid_medianint_YFP = YFP_MedianIntensity_nuc_tot{cell,1}(valid_indices_YFP, 2);

    %% bg rescale
    % % %     valid_medianint_YFP = YFP_MedianIntensity_nuc_tot{cell,1}(valid_indices_YFP, 2);
    % %     valid_medianint_YFP_bg_rescale = (valid_medianint_YFP-YFP_bg)./250; %% set bg to 0, every 1000 increase scale to 1
    % % %     plot(YFP_MedianIntensity_nuc_tot{cell,1}(valid_indices_YFP, 1), valid_medianint_YFP_bg_rescale, 'LineWidth', 2);

    valid_medianint_YFP_bg_rescale = (valid_medianint_YFP-YFP_bg); %% set bg to 0


    valid_time_int_YFP = [valid_time_YFP valid_medianint_YFP_bg_rescale];
    order = 3;
    framelen = 11;

    YFP_sgf = sgolayfilt(valid_time_int_YFP,order,framelen);

    YFP_sgf_tot{tot_cell_num} = YFP_sgf;

    % %      [minValue, minIndex] = min(valid_time_YFP);
    % %     YFP_basal_amp(tot_cell_num,1) = valid_medianint_YFP_bg_rescale(minIndex, 1);

%     [minValue, minIndex] = min(valid_time_YFP);
% 
%     YFP_basal_time(tot_cell_num,1)  = valid_time_YFP(minIndex);
%     YFP_basal_amp(tot_cell_num,1)  =YFP_sgf(minIndex, 2);

        [minValue, minIndex] = min(valid_time_YFP(valid_time_YFP>=21));
         YFP_basal_time(tot_cell_num,1)  = valid_time_YFP(valid_time_YFP==minValue);
        YFP_basal_amp(tot_cell_num,1)  =YFP_sgf(valid_time_YFP==minValue, 2);

    %% visualization
    plot(valid_time_int_YFP(:,1),valid_time_int_YFP(:,2),':','LineWidth',2)
    hold on
    plot(YFP_sgf(:,1),YFP_sgf(:,2),'.-','LineWidth',2)
    hold on

    tot_cell_num = tot_cell_num +1;
end

YFP_peak_amp = zeros(tot_cell_num-1,1);
YFP_peak_time = zeros(tot_cell_num-1,1);
%% extract the peak(s)
for c_num = 1: (tot_cell_num-1)

    YFP_sgf_time = YFP_sgf_tot{c_num}(:,1);
    YFP_sgf_int = YFP_sgf_tot{c_num}(:,2);
    [YFP_pks,locs] = findpeaks(YFP_sgf_int(YFP_sgf_time<602));

    [YFP_max_peak, YFP_max_idx] = max(YFP_pks); % Find the highest peak
    YFP_peak_location = locs(YFP_max_idx);

    YFP_peak_amp_ori(c_num,1) = YFP_sgf_tot{c_num}(YFP_peak_location,2);
    YFP_peak_amp_thresh(c_num,1)=0.9.*YFP_peak_amp_ori(c_num,1);
   later_basal = find(YFP_sgf_tot{c_num}(:,1)>=YFP_basal_time(c_num,1));
    idx_peak = find(YFP_sgf_tot{c_num}(later_basal,2) >= YFP_peak_amp_thresh(c_num,1), 1, 'first');
    YFP_peak_time(c_num,1)  = YFP_sgf_tot{c_num}(later_basal(idx_peak),1);
    YFP_peak_amp(c_num,1)  = YFP_sgf_tot{c_num}(later_basal(idx_peak),2);
% %     
% %     YFP_peak_amp(c_num,1) = YFP_sgf_tot{c_num}(YFP_peak_location,2);
% %     YFP_peak_time(c_num,1) = YFP_sgf_tot{c_num}(YFP_peak_location,1);
end

save('/Users/yijiachen/Documents/B_Cell_project/HVN029/data_saved/HVN029_YFP_nuc_peak_amp_time_90_perc.mat','YFP_peak_amp','YFP_peak_time','YFP_basal_amp','YFP_basal_time')


%% visualization
scatter(YFP_peak_time,YFP_peak_amp,100,'filled')
scatter(YFP_basal_time,YFP_basal_amp,100,'filled')


YFP_rgb_colors = [255,204,0; 255,255,0; 204,153,0; 255,153,0; 255,255,102; 184,180,0; 125,122,0; 204,204,0; 255 255 51];
YFP_colors = YFP_rgb_colors / 255;
box off; % Turn off the box surrounding the plot
ax = gca;
ax.ColorOrder = YFP_colors;
ax.TickDir = 'in'; % Set the direction of the ticks

% ylim([650 950])
% ylim([0 1])
xlim([0 960])
xticks([1 60:60:960])
xticklabels({'0', '1', '2', '3', '4', '5', '6','7', '8', '9', '10', '11', '12', '13','14', '15', '16'});
title('YFP nuclear median intensity')
xlabel('time (h)')
ylabel('YFP nuclear median intensity')
legend('signal','sgolay')
% saveas(YFP_nuclear_medianint_plot_smoothed,'/Users/yijiachen/Documents/B_Cell_project/HVN029/images/YFP_nuc_rep_traj/fluor_plot_medianint_smoothed_90_perc.png')
% saveas(YFP_nuclear_medianint_plot_smoothed,sprintf('/Users/yijiachen/Documents/B_Cell_project/HVN029/images/YFP_nuc_rep_traj/fluor_plot_medianint_smoothed_c%d.png',cell))




