for s = 1:9
    time_point=1;
    for t = [1 2 3 4 7 13 19 25 49 97 145]

        path = sprintf('/Users/yijiachen/Documents/B_Cell_project/panel_cd/HVN032_live_nucleus_mask/CFP_masked/s%d', s);
        addpath(path);
        filename_CFP_nlm = ['CFP_img_masked_', num2str(t), '.tiff'];
        info = imfinfo(filename_CFP_nlm);

        CFP_img_nlm = zeros(info(1).Height,info(1).Width,length(info),'uint16');
        for read = 1:length(info)
            CFP_img_nlm(:,:,read) = imread(filename_CFP_nlm,read);
        end
        label_CFP_nlm = bwlabel(CFP_img_nlm);

        CFP_stats = regionprops("table",label_CFP_nlm,CFP_img_nlm,"Centroid",'MeanIntensity',"MajorAxisLength","MinorAxisLength") ;

        diameters = mean([CFP_stats.MajorAxisLength CFP_stats.MinorAxisLength],2);
        radii_CFP = diameters/2;

        tabledata_CFP = table2array(CFP_stats);
        CFP_nucleus_live_mask_all_t.(['s', num2str(s)]).('nlm')(1:size(CFP_stats,1),time_point)=tabledata_CFP(:,5);
        CFP_nucleus_live_mask_all_t.(['s', num2str(s)]).('nlm_centroid')(1:size(CFP_stats,1),((2*(time_point-1)+1):(2*(time_point-1)+1)+1))=tabledata_CFP(:,1:2);
        CFP_nucleus_live_mask_all_t.(['s', num2str(s)]).('radius')(1:size(CFP_stats,1),time_point)=radii_CFP;
%         %% visualization
%         centroid_y = CFP_nucleus_live_mask_all_t.(['s', num2str(s)]).('nlm_centroid')(1:size(CFP_stats,1),((2*(time_point-1)+1)));
%         centroid_x = CFP_nucleus_live_mask_all_t.(['s', num2str(s)]).('nlm_centroid')(1:size(CFP_stats,1),((2*(time_point-1)+1)+1));
%         figure
%         imagesc(CFP_img_nlm)
%         truesize([500 500]);
%         colormap gray
%         colorbar
%         title("CFP_img_mask")
%         hold on
%         scatter(centroid_y, centroid_x, 'r.', 'SizeData', 100); % Adjust 'SizeData' as needed
% 



        path = sprintf('/Users/yijiachen/Documents/B_Cell_project/panel_cd/HVN032_live_nucleus_mask/YFP_masked/s%d', s);
        addpath(path);
        filename_YFP_nlm = ['YFP_img_masked_', num2str(t), '.tiff'];
        info = imfinfo(filename_YFP_nlm);

        YFP_img_nlm = zeros(info(1).Height,info(1).Width,length(info),'uint16');
        for read = 1:length(info)
            YFP_img_nlm(:,:,read) = imread(filename_YFP_nlm,read);
        end
        label_YFP_nlm = bwlabel(YFP_img_nlm);

        YFP_stats = regionprops("table",label_YFP_nlm,YFP_img_nlm,"Centroid",'MeanIntensity',"MajorAxisLength","MinorAxisLength") ;
        centroid_nucl = YFP_stats.Centroid;
        diameters = mean([YFP_stats.MajorAxisLength YFP_stats.MinorAxisLength],2);
        radii_YFP = diameters/2;

        tabledata_YFP = table2array(YFP_stats);
        YFP_nucleus_live_mask_all_t.(['s', num2str(s)]).('nlm')(1:size(YFP_stats,1),time_point)=tabledata_YFP(:,5);
        YFP_nucleus_live_mask_all_t.(['s', num2str(s)]).('nlm_centroid')(1:size(YFP_stats,1),((2*(time_point-1)+1):(2*(time_point-1)+1)+1))=tabledata_YFP(:,1:2);
        YFP_nucleus_live_mask_all_t.(['s', num2str(s)]).('radius')(1:size(YFP_stats,1),time_point)=radii_YFP;

%         close all
        time_point=time_point+1;

    end
    save('/Users/yijiachen/Documents/B_Cell_project/panel_cd/data_saved/HVN032_CFP_nucleus_live_mask_value.mat','CFP_nucleus_live_mask_all_t')
    save('/Users/yijiachen/Documents/B_Cell_project/panel_cd/data_saved/HVN032_YFP_nucleus_live_mask_value.mat','YFP_nucleus_live_mask_all_t')

end