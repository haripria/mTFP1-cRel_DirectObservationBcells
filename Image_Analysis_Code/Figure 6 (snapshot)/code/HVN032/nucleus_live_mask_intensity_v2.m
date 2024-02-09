for s = 1:9
    time_point=1;
    for t =[1 2 3 4 7 13 19 25 49 97 145]
        %% load live nucleus mask
        path = sprintf('/Users/yijiachen/Documents/B_Cell_project/panel_cd/HVN032_live_nucleus_mask/live_nucleus_mask/s%d', s);
        addpath(path);
        filename_final_live_mask = ['live_nucleus_mask_', num2str(t,'%01.f'), '.tiff'];
        final_live_mask=imread(filename_final_live_mask);

        %% load CFP image
        path_2 = sprintf('/Users/yijiachen/Documents/B_Cell_project/panel_cd/HVN032_Fluo_9_scenes/HVN032_Fluo_ALL_s%d', s);
        addpath(path_2);
        filename_CFP = ['HVN032_Fluo_ALL_t', num2str(t,'%03.f'), 'c5_ORG', '.tif'];
        info = imfinfo(filename_CFP);
        length(info);
        info(1).Width;
        info(1).Height;
        info(1).BitDepth;
        CFP_img = zeros(info(1).Height,info(1).Width,length(info),'uint16');
        for i = 1:length(info)
            CFP_img(:,:,i) = imread(filename_CFP,i);
        end

        CFP_img_nlm = CFP_img;
        CFP_img_nlm(~final_live_mask) = 0;

        %         Save the image as a TIFF file without reducing the resolution
        savePath = sprintf('/Users/yijiachen/Documents/B_Cell_project/panel_cd/HVN032_live_nucleus_mask/CFP_masked/s%d/CFP_img_masked_%d.tiff', s, t);
        imwrite(CFP_img_nlm, savePath);


        CFP_stats = regionprops("table",final_live_mask,CFP_img_nlm,'Area',"Centroid",'MeanIntensity',"MajorAxisLength","MinorAxisLength",'PixelValues') ;
        CFP_stats = CFP_stats(CFP_stats.Area>=100 & CFP_stats.Area <= 800,:);
        diameters = mean([CFP_stats.MajorAxisLength CFP_stats.MinorAxisLength],2);
        radii_CFP = diameters/2;
        fn=['s', num2str(s)];

        for row_num = 1:size(CFP_stats,1)
            CFP_nucleus_live_mask_all_t.(fn).('nlm_medianint')(row_num,time_point)=median(cell2mat(CFP_stats.PixelValues(row_num,1)));
            CFP_nucleus_live_mask_all_t.(fn).('nlm_sumint')(row_num,time_point)=sum(cell2mat(CFP_stats.PixelValues(row_num,1)));
            CFP_nucleus_live_mask_all_t.(fn).('nlm_meanint')(row_num,time_point)=CFP_stats.MeanIntensity(row_num,1);
            CFP_nucleus_live_mask_all_t.(fn).('nlm_centroid')(row_num,((2*(time_point-1)+1):(2*(time_point-1)+1)+1))=CFP_stats.Centroid(row_num,1:2);
            CFP_nucleus_live_mask_all_t.(fn).('radius')(row_num,time_point)=radii_CFP(row_num,1);
            CFP_nucleus_live_mask_all_t.(fn).('area')(row_num,time_point)=CFP_stats.Area(row_num,1);
        end

      
        %% load YFP image
        filename_YFP = ['HVN032_Fluo_ALL_t', num2str(t,'%03.f'), 'c4_ORG', '.tif'];
        info = imfinfo(filename_YFP);
        length(info);
        info(1).Width;
        info(1).Height;
        info(1).BitDepth;
        YFP_img = zeros(info(1).Height,info(1).Width,length(info),'uint16');
        for i = 1:length(info)
            YFP_img(:,:,i) = imread(filename_YFP,i);
        end

        YFP_img_nlm = YFP_img;
        YFP_img_nlm(~final_live_mask) = 0;

        % Save the image as a TIFF file without reducing the resolution
        savePath = sprintf('/Users/yijiachen/Documents/B_Cell_project/panel_cd/HVN032_live_nucleus_mask/YFP_masked/s%d/YFP_img_masked_%d.tiff', s, t);
        imwrite(YFP_img_nlm, savePath);


        YFP_stats = regionprops("table",final_live_mask,YFP_img_nlm,'Area',"Centroid",'MeanIntensity',"MajorAxisLength","MinorAxisLength",'PixelValues') ;
        YFP_stats = YFP_stats(YFP_stats.Area>=100 & YFP_stats.Area <= 800,:);
        diameters = mean([YFP_stats.MajorAxisLength YFP_stats.MinorAxisLength],2);
        radii_YFP = diameters/2;

        for row_num_2 = 1:size(YFP_stats,1)
            YFP_nucleus_live_mask_all_t.(fn).('nlm_medianint')(row_num_2,time_point)=median(cell2mat(YFP_stats.PixelValues(row_num_2,1)));
            YFP_nucleus_live_mask_all_t.(fn).('nlm_sumint')(row_num_2,time_point)=sum(cell2mat(YFP_stats.PixelValues(row_num_2,1)));
            YFP_nucleus_live_mask_all_t.(fn).('nlm_meanint')(row_num_2,time_point)=YFP_stats.MeanIntensity(row_num_2,1);
            YFP_nucleus_live_mask_all_t.(fn).('nlm_centroid')(row_num_2,((2*(time_point-1)+1):(2*(time_point-1)+1)+1))=YFP_stats.Centroid(row_num_2,1:2);
            YFP_nucleus_live_mask_all_t.(fn).('radius')(row_num_2,time_point)=radii_YFP(row_num_2,1);
            YFP_nucleus_live_mask_all_t.(fn).('area')(row_num,time_point)=YFP_stats.Area(row_num,1);
        end

        time_point=time_point+1;


%           %% visualization
%                 figure
%                 imagesc(CFP_img_nlm)
%                 truesize([500 500]);
%                 colormap gray
%                 colorbar
%                 title("CFP_img_nlm")
%                 hold on
%                 scatter(CFP_stats.Centroid(:, 1), CFP_stats.Centroid(:, 2), 'r.', 'SizeData', 100); % Adjust 'SizeData' as needed
%                 % Plot the mean intensity values
%                 outlier_count=1;
%                 for i = 1:size(CFP_stats,1)
%                     x = CFP_stats.Centroid(i, 1); % X-coordinate of the centroid
%                     y = CFP_stats.Centroid(i, 2); % Y-coordinate of the centroid
% %                     meanIntensity = floor(CFP_stats.MeanIntensity(i)); % Mean intensity value
% %                     % Add the mean intensity value as text at the centroid
% %                     text(x, y+15, num2str(meanIntensity), 'Color', 'yellow', 'FontSize', 10, 'FontWeight', 'bold');
% %                     areavalue = floor(CFP_stats.Area(i)); % Mean intensity value
% %                     % Add the mean intensity value as text at the centroid
% %                     text(x, y+15, num2str(areavalue), 'Color', 'yellow', 'FontSize', 10, 'FontWeight', 'bold');
%                     pixelvalue_CFP = median(cell2mat(CFP_stats.PixelValues(i))); % Mean intensity value
% 
%                     
%                     if pixelvalue_CFP >=3500
%                     % Add the mean intensity value as text at the centroid
%                     text(x, y+15, num2str(pixelvalue_CFP), 'Color', 'yellow', 'FontSize', 10, 'FontWeight', 'bold');
%                     outlier_count=outlier_count+1;
%                     end
%                     
%                 end
%                 outlier_count
% %                 figure
% %                 scatter(1:size(pixelvalue,1),pixelvalue)
% 
% 
%                 figure
%                 imagesc(YFP_img_nlm)
%                 truesize([500 500]);
%                 colormap gray
%                 colorbar
%                 title("YFP_img_nlm")
%                 hold on
%                 scatter(YFP_stats.Centroid(:, 1), YFP_stats.Centroid(:, 2), 'r.', 'SizeData', 100); % Adjust 'SizeData' as needed
%                 % Plot the mean intensity values
%                 for i = 1:size(YFP_stats,1)
%                     x = YFP_stats.Centroid(i, 1); % X-coordinate of the centroid
%                     y = YFP_stats.Centroid(i, 2); % Y-coordinate of the centroid
% %                     meanIntensity = floor(CFP_stats.MeanIntensity(i)); % Mean intensity value
% %                     % Add the mean intensity value as text at the centroid
% %                     text(x, y+15, num2str(meanIntensity), 'Color', 'yellow', 'FontSize', 10, 'FontWeight', 'bold');
% %                     areavalue = floor(CFP_stats.Area(i)); % Mean intensity value
% %                     % Add the mean intensity value as text at the centroid
% %                     text(x, y+15, num2str(areavalue), 'Color', 'yellow', 'FontSize', 10, 'FontWeight', 'bold');
%                     
%                     pixelvalue_YFP = median(cell2mat(YFP_stats.PixelValues(i))); % Mean intensity value
%                     % Add the mean intensity value as text at the centroid
%                     if pixelvalue_YFP >=1000
%                     text(x, y+15, num2str(pixelvalue_YFP), 'Color', 'yellow', 'FontSize', 10, 'FontWeight', 'bold');
%                     end
%                 end

    end
    save('/Users/yijiachen/Documents/B_Cell_project/panel_cd/data_saved/HVN032_CFP_nucleus_live_mask_value.mat','CFP_nucleus_live_mask_all_t')
    save('/Users/yijiachen/Documents/B_Cell_project/panel_cd/data_saved/HVN032_YFP_nucleus_live_mask_value.mat','YFP_nucleus_live_mask_all_t')

end