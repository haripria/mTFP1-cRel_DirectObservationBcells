
load('/Users/yijiachen/Documents/B_Cell_project/HVN037/data_saved/HVN037_CFP_YFP_nuc_median_mean_int.mat')
load('/Users/yijiachen/Documents/B_Cell_project/HVN037/data_saved/HVN037_CFP_YFP_wc_median_mean_int.mat')

%% nuc
tot_cell_num = size(CFP_MedianIntensity_nuc_tot,1);
CFP_nuc_mn = zeros(tot_cell_num,1);
CFP_nuc_sd = zeros(tot_cell_num,1);

YFP_nuc_mn = zeros(tot_cell_num,1);
YFP_nuc_sd = zeros(tot_cell_num,1);


for cn = 1:tot_cell_num
    CFP_nuc_columnData = CFP_MedianIntensity_nuc_tot{cn,1}(1:37,2);
    CFP_nuc_mn(cn,1)=mean(CFP_nuc_columnData(~isnan(CFP_nuc_columnData)));
    CFP_nuc_sd(cn,1)=std(CFP_nuc_columnData(~isnan(CFP_nuc_columnData)));


    YFP_nuc_columnData = YFP_MedianIntensity_nuc_tot{cn,1}(1:37,2);
    YFP_nuc_mn(cn,1)=mean(YFP_nuc_columnData(~isnan(YFP_nuc_columnData)));
    YFP_nuc_sd(cn,1)=std(YFP_nuc_columnData(~isnan(YFP_nuc_columnData)));

end
CFP_nuc_cv=CFP_nuc_sd./CFP_nuc_mn;
YFP_nuc_cv=YFP_nuc_sd./YFP_nuc_mn;

save('/Users/yijiachen/Documents/B_Cell_project/HVN037/data_saved/HVN037_CFP_YFP_nuc_cv.mat','CFP_nuc_cv','YFP_nuc_cv')

%% wc
tot_cell_num = size(CFP_ConcIntensity_wc_tot,1);
CFP_wc_mn = zeros(tot_cell_num,1);
CFP_wc_sd = zeros(tot_cell_num,1);

YFP_wc_mn = zeros(tot_cell_num,1);
YFP_wc_sd = zeros(tot_cell_num,1);


for cn = 1:tot_cell_num
    CFP_wc_columnData = CFP_ConcIntensity_wc_tot{cn,1}(1:37,2);
    CFP_wc_mn(cn,1)=mean(CFP_wc_columnData(~isnan(CFP_wc_columnData)));
    CFP_wc_sd(cn,1)=std(CFP_wc_columnData(~isnan(CFP_wc_columnData)));


    YFP_wc_columnData = YFP_ConcIntensity_wc_tot{cn,1}(1:37,2);
    YFP_wc_mn(cn,1)=mean(YFP_wc_columnData(~isnan(YFP_wc_columnData)));
    YFP_wc_sd(cn,1)=std(YFP_wc_columnData(~isnan(YFP_wc_columnData)));

end
CFP_wc_cv=CFP_wc_sd./CFP_wc_mn;
YFP_wc_cv=YFP_wc_sd./YFP_wc_mn;

save('/Users/yijiachen/Documents/B_Cell_project/HVN037/data_saved/HVN037_CFP_YFP_wc_cv.mat','CFP_wc_cv','YFP_wc_cv')


%% CFP nuc vs wc CV Violin Plot

CFP_nuc_data = [repmat({'Nuc'}, length(CFP_nuc_cv), 1)];
CFP_wc_data = [repmat({'WC'}, length(CFP_wc_cv), 1)];


CFP_CV_violin = figure;

swarmchart(categorical(CFP_nuc_data), CFP_nuc_cv, 30,'filled','MarkerFaceColor','#00CCFF')
hold on
swarmchart(categorical(CFP_wc_data),CFP_wc_cv,30,'filled','MarkerFaceColor','#0000FF')

xlabel('X-axis Label', 'FontSize', 12);
ylabel('Y-axis Label', 'FontSize', 12);
title('CFP CV Violin Plot', 'FontSize', 14);
ax = gca;
ax.FontSize = 10;
ax.XGrid = 'on';
ax.YGrid = 'on';
set(CFP_CV_violin, 'Position', [0, 0, 123, 370]); 

saveas(CFP_CV_violin,'/Users/yijiachen/Documents/B_Cell_project/HVN037/images/CV_violin/CFP_CV_medianint_violin.png')


%% YFP nuc vs wc CV Violin Plot

YFP_nuc_data = [repmat({'Nuc'}, length(YFP_nuc_cv), 1)];
YFP_wc_data = [repmat({'WC'}, length(YFP_wc_cv), 1)];


YFP_CV_violin = figure;

swarmchart(categorical(YFP_nuc_data), YFP_nuc_cv, 30, 'filled','MarkerFaceColor','#FFCC00');
hold on
swarmchart(categorical(YFP_wc_data),YFP_wc_cv,30,'filled','MarkerFaceColor','#B8B400')

xlabel('X-axis Label', 'FontSize', 12);
ylabel('Y-axis Label', 'FontSize', 12);
title('YFP CV Violin Plot', 'FontSize', 14);
ax = gca;
ax.FontSize = 10;
ax.XGrid = 'on';
ax.YGrid = 'on';
set(YFP_CV_violin, 'Position', [0, 0, 123, 370]); 

saveas(YFP_CV_violin,'/Users/yijiachen/Documents/B_Cell_project/HVN037/images/CV_violin/YFP_CV_medianint_violin.png')

% % %% CFP YFP CV Violin Plot
% %
% % CFP_x_data = [repmat({'CFP'}, length(CFP_nuc_cv), 1)];
% % YFP_x_data = [repmat({'YFP'}, length(YFP_nuc_cv), 1)];
% %
% %
% % figure;
% %
% % swarmchart(categorical(CFP_x_data), CFP_nuc_cv, 50, 'filled', 'blue');
% % hold on
% % swarmchart(categorical(YFP_x_data),YFP_nuc_cv,50,'filled','yellow')
% %
% % xlabel('X-axis Label', 'FontSize', 12);
% % ylabel('Y-axis Label', 'FontSize', 12);
% % title('CV Violin Plot', 'FontSize', 14);
% % ax = gca;
% % ax.FontSize = 10;
% % ax.XGrid = 'on';
% % ax.YGrid = 'on';
% %
% %
% % %% CFP YFP wc CV Violin Plot
% %
% % CFP_x_data = [repmat({'CFP'}, length(CFP_wc_cv), 1)];
% % YFP_x_data = [repmat({'YFP'}, length(YFP_wc_cv), 1)];
% %
% %
% % figure;
% %
% % swarmchart(categorical(CFP_x_data), CFP_wc_cv, 50, 'filled', 'blue');
% % hold on
% % swarmchart(categorical(YFP_x_data),YFP_wc_cv,50,'filled','yellow')
% %
% % xlabel('X-axis Label', 'FontSize', 12);
% % ylabel('Y-axis Label', 'FontSize', 12);
% % title('CV Violin Plot', 'FontSize', 14);
% % ax = gca;
% % ax.FontSize = 10;
% % ax.XGrid = 'on';
% % ax.YGrid = 'on';

% %% CFP YFP CV histogram
%
% figure;
% histogram(CFP_nuc_cv, 'NumBins', 20, 'FaceColor', [0, 0.5, 1]);
% xlabel('CFP CV Value', 'FontSize', 12);
% ylabel('Frequency', 'FontSize', 12);
% title('Histogram of CFP Coefficient of Variation', 'FontSize', 14);
% ylim([0,25])
% ax = gca;
% ax.FontSize = 10;
% ax.XGrid = 'on';
% ax.YGrid = 'on';
%
% figure;
% histogram(YFP_nuc_cv, 'NumBins', 20, 'FaceColor', [0, 0.5, 1]);
% xlabel('YFP CV Value', 'FontSize', 12);
% ylabel('Frequency', 'FontSize', 12);
% title('Histogram of YFP Coefficient of Variation', 'FontSize', 14);
% ylim([0,25])
% ax = gca;
% ax.FontSize = 10;
% ax.XGrid = 'on';
% ax.YGrid = 'on';
