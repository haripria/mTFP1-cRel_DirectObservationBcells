load('/Users/yijiachen/Documents/B_Cell_project/HVN029/data_saved/HVN029_CFP_nuc_peak_amp_time.mat')
CFP_pb_diff_wt = CFP_pb_diff;
CFP_fold_change_wt = CFP_fold_change;
clear CFP_pb_diff CFP_fold_change


load('/Users/yijiachen/Documents/B_Cell_project/HVN032/data_saved/HVN032_CFP_nuc_peak_amp_time.mat')
CFP_pb_diff_eko = CFP_pb_diff;
CFP_fold_change_eko = CFP_fold_change;
clear CFP_pb_diff CFP_fold_change



pb_diff_hist = figure;
histogram(CFP_pb_diff_wt,10,'FaceColor','black','FaceAlpha',0.5)
hold on 
histogram(CFP_pb_diff_eko,10,'FaceColor','#FF0066','FaceAlpha',0.5)
title('cRel peak basal diff histogram')
legend('WT','eKO')

saveas(pb_diff_hist,'/Users/yijiachen/Documents/B_Cell_project/images/CFP_pb_diff_wt_eKO_histogram.png')

fold_change_hist = figure;
histogram(CFP_fold_change_wt,10,'FaceColor','black','FaceAlpha',0.5)
hold on
histogram(CFP_fold_change_eko,10,'FaceColor','#FF0066','FaceAlpha',0.5)
title('cRel fold change histogram')
legend('WT','eKO')

saveas(fold_change_hist,'/Users/yijiachen/Documents/B_Cell_project/images/CFP_fold_change_wt_eKO_histogram.png')

% load('/Users/yijiachen/Documents/B_Cell_project/HVN029/data_saved/HVN029_YFP_nuc_peak_amp_time.mat')
% load('/Users/yijiachen/Documents/B_Cell_project/HVN032/data_saved/HVN032_YFP_nuc_peak_amp_time.mat')