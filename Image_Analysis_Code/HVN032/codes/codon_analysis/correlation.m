load ('/Users/yijiachen/Documents/B_Cell_project/HVN032/data_saved/HVN032_CFP_nuc_peak_amp_time.mat')
load ('/Users/yijiachen/Documents/B_Cell_project/HVN032/data_saved/HVN032_YFP_nuc_peak_amp_time.mat')

% 'CFP_peak_amp','CFP_peak_time','CFP_basal_amp','CFP_fold_change'

[r,pval] = corr(CFP_fold_change,YFP_peak_amp)
r_2 = r^2

[r,pval] = corr(YFP_peak_amp,CFP_fold_change)
r_2 = r^2


[r,pval] = corr(CFP_peak_amp,CFP_peak_time)
r_2 = r^2

[r,pval] = corr(CFP_peak_amp,CFP_basal_amp)
r_2 = r^2

[r,pval] = corr(CFP_peak_amp,CFP_fold_change)
r_2 = r^2

[r,pval] = corr(CFP_peak_time,CFP_basal_amp)
r_2 = r^2

[r,pval] = corr(CFP_peak_time,CFP_fold_change)
r_2 = r^2

[r,pval] = corr(CFP_basal_amp,CFP_fold_change)
r_2 = r^2


[r,pval] = corr(YFP_peak_amp,YFP_peak_time)
r_2 = r^2

[r,pval] = corr(YFP_peak_amp,YFP_basal_amp)
r_2 = r^2

% [r,pval] = corr(YFP_peak_amp,YFP_fold_change)
% r_2 = r^2

[r,pval] = corr(YFP_peak_time,YFP_basal_amp)
r_2 = r^2

% [r,pval] = corr(YFP_peak_time,YFP_fold_change)
% r_2 = r^2

% [r,pval] = corr(YFP_basal_amp,YFP_fold_change)
% r_2 = r^2

