load('/Users/yijiachen/Documents/B_Cell_project/HVN032/data_saved/HVN032_CFP_YFP_nuc_median_mean_int.mat')


%% CFP nuc median int traj
CFP_bg = 2200;

CFP_nuclear_medianint_plot = figure;
s1_num = 15;
s4_num = 8;
s5_num = 5;
s6_num = 4;
s7_num = 7;
s8_num = 8;
s9_num = 6;

s1_cell = [1,2,3,4,5,6,7,8,9,15];
s4_cell = [1,3,5,6];
s5_cell = [1,2,3];
s6_cell = [1,3,4];
s7_cell = [2,3,4,5];
s8_cell = [4,5,6];
s9_cell = [1,3];

tot_cell = [s1_cell, s4_cell+s1_num, s5_cell+s1_num+s4_num, s6_cell+s1_num+s4_num+s5_num, s7_cell+s1_num+s4_num+s5_num+s6_num, s8_cell+s1_num+s4_num+s5_num+s6_num+s7_num, s9_cell+s1_num+s4_num+s5_num+s6_num+s7_num+s8_num];

% for s = 1

for cell = tot_cell
    %[1,2,3,4,6,7,8,9,15]
    %% CFP nuclear median intensity traj

    valid_indices_CFP = ~isnan(CFP_MedianIntensity_nuc_tot{cell,1}(:, 2)) & CFP_MedianIntensity_nuc_tot{cell,1}(:, 2) ~= 0;
%% bg rescale
    valid_medianint_CFP = CFP_MedianIntensity_nuc_tot{cell,1}(valid_indices_CFP, 2);
    valid_medianint_CFP_bg_rescale = (valid_medianint_CFP-CFP_bg)./1000; %% set bg to 0, every 1000 increase scale to 1
    plot(CFP_MedianIntensity_nuc_tot{cell,1}(valid_indices_CFP, 1), valid_medianint_CFP_bg_rescale, 'LineWidth', 2);

% %     %     CFP_MedianIntensity_nuc_tot_bg_rescaled = (CFP_MedianIntensity_nuc_tot{cell,1}(valid_indices_CFP, 2)- CFP_bg_rescale_lowb)/(CFP_bg_rescale_upb-CFP_bg_rescale_lowb)*(new_max-new_min);
% %     plot(CFP_MedianIntensity_nuc_tot{cell,1}(valid_indices_CFP, 1), CFP_MedianIntensity_nuc_tot{cell,1}(valid_indices_CFP, 2), 'LineWidth', 2);
% %     %     plot(CFP_MedianIntensity_nuc_tot{cell,1}(valid_indices_CFP, 1), CFP_MedianIntensity_nuc_tot_bg_rescaled, 'LineWidth', 2);

    hold on

end

CFP_rgb_colors = [0,204,255; 0,0,255; 0,102,255; 0,153,255; 153,204,255; 51,51,204; 51,102,204; 0,255,255; 102,255,255];
CFP_colors = CFP_rgb_colors / 255;

box off; % Turn off the box surrounding the plot
ax = gca;
ax.ColorOrder = CFP_colors;
ax.TickDir = 'in'; % Set the direction of the ticks

ylim([0 1.5])
% ylim([2000 3200])
xlim([0 960])
xticks([1 60:60:960])
xticklabels({'0', '1', '2', '3', '4', '5', '6','7', '8', '9', '10', '11', '12', '13','14', '15', '16'});
title('CFP nuclear median intensity')
xlabel('time (h)')
ylabel('CFP nuclear median intensity')
% end
saveas(CFP_nuclear_medianint_plot,'/Users/yijiachen/Documents/B_Cell_project/HVN032/images/CFP_nuc_rep_traj/fluor_plot_medianint.png')


%% YFP nuc median int traj
YFP_bg = 650;

YFP_nuclear_medianint_plot = figure;

for cell = tot_cell
    %% YFP nuclear median intensity traj

    valid_indices_YFP = ~isnan(YFP_MedianIntensity_nuc_tot{cell,1}(:, 2)) & YFP_MedianIntensity_nuc_tot{cell,1}(:, 2) ~= 0;
    
     %% bg rescale
    valid_medianint_YFP = YFP_MedianIntensity_nuc_tot{cell,1}(valid_indices_YFP, 2);
    valid_medianint_YFP_bg_rescale = (valid_medianint_YFP - YFP_bg)./250; %% set bg to 0, every 250 increase scale to 1
    plot(YFP_MedianIntensity_nuc_tot{cell,1}(valid_indices_YFP, 1), valid_medianint_YFP_bg_rescale, 'LineWidth', 2);

% %     plot(YFP_MedianIntensity_nuc_tot{cell,1}(valid_indices_YFP, 1), YFP_MedianIntensity_nuc_tot{cell,1}(valid_indices_YFP, 2), 'LineWidth', 2);
% % 
% %     %     YFP_MedianIntensity_nuc_tot_bg_rescaled = (YFP_MedianIntensity_nuc_tot{cell,1}(valid_indices_YFP, 2)- YFP_bg_rescale_lowb)/(YFP_bg_rescale_upb-YFP_bg_rescale_lowb)*(new_max-new_min);
% %     %     %         plot(YFP_MedianIntensity_nuc_tot{cell,1}(valid_indices_YFP, 1), YFP_MedianIntensity_nuc_tot{cell,1}(valid_indices_YFP, 2), 'LineWidth', 2);
% %     %     plot(YFP_MedianIntensity_nuc_tot{cell,1}(valid_indices_YFP, 1), YFP_MedianIntensity_nuc_tot_bg_rescaled, 'LineWidth', 2);

    hold on

end

YFP_rgb_colors = [255,204,0; 255,255,0; 204,153,0; 255,153,0; 255,255,102; 184,180,0; 125,122,0; 204,204,0; 255 255 51];
YFP_colors = YFP_rgb_colors / 255;


box off; % Turn off the box surrounding the plot
ax = gca;
ax.ColorOrder = YFP_colors;
ax.TickDir = 'in'; % Set the direction of the ticks


ylim([0 1.5])
% ylim([650 950])
xlim([0 960])
xticks([1 60:60:960])
xticklabels({'0', '1', '2', '3', '4', '5', '6','7', '8', '9', '10', '11', '12', '13','14', '15', '16'});

title('YFP nuclear median intensity')
xlabel('time (h)')
ylabel('YFP nuclear median intensity')
% legend ('YFP','YFP')
% end

saveas(YFP_nuclear_medianint_plot,'/Users/yijiachen/Documents/B_Cell_project/HVN032/images/YFP_nuc_rep_traj/fluor_plot_medianint.png')

save('/Users/yijiachen/Documents/B_Cell_project/HVN032/data_saved/HVN032_s1_s4_to_s9_tot_cell.mat','tot_cell')
