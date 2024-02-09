load('/Users/yijiachen/Documents/B_Cell_project/panel_cd/data_saved/HVN029_CFP_bg_mean_median.mat')
load('/Users/yijiachen/Documents/B_Cell_project/panel_cd/data_saved/HVN032_CFP_bg_mean_median.mat')

CFP_bg_mean = (bg_mean_HVN029 + bg_mean_HVN032)/2;
CFP_bg_median = (bg_median_HVN029 + bg_median_HVN032)/2;
save('/Users/yijiachen/Documents/B_Cell_project/panel_cd/data_saved/CFP_bg_mean_median.mat','CFP_bg_mean','CFP_bg_median')
clear

load('/Users/yijiachen/Documents/B_Cell_project/panel_cd/data_saved/HVN029_YFP_bg_mean_median.mat')
load('/Users/yijiachen/Documents/B_Cell_project/panel_cd/data_saved/HVN032_YFP_bg_mean_median.mat')

YFP_bg_mean = (bg_mean_HVN029 + bg_mean_HVN032)/2;
YFP_bg_median = (bg_median_HVN029 + bg_median_HVN032)/2;
save('/Users/yijiachen/Documents/B_Cell_project/panel_cd/data_saved/YFP_bg_mean_median.mat','YFP_bg_mean','YFP_bg_median')
