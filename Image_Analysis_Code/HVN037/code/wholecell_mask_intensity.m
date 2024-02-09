for s = [2,5,6,7,8]
%     [1,2,5,6,7,8]
    time_point=1;
    for t = 1:92
        %1:150
        %[1 2 3 4 7 13 19 25 49 97 145]

        t_wc = (t-1)*10+1;
        %% load wholecell mask
        path = sprintf('/Volumes/Data/Mark/HVN037a_Masks/s%d', s);
        addpath(path);

        filename_wc = sprintf('HVN037a_NoStim_TFP_Venus_Cherry_IkBeKO_40XOil_RedBF-02_scene%d_t%03d_ORG_cp_masks.png', s, t_wc);
        wc_mask=imread(filename_wc);

        %% load CFP image
        path_2 = sprintf('/Volumes/Data/Mark/HVN037a_NoStim_TFP_Venus_Cherry_IkBeKO_40XOil_RedBF-02_scene%d', s);
        addpath(path_2);

        filename_CFP = sprintf('HVN037a_NoStim_TFP_Venus_Cherry_IkBeKO_40XOil_RedBF-02_scene%d_t%02dCFP_3_ORG.tif', s, t);

        info = imfinfo(filename_CFP);
        length(info);
        info(1).Width;
        info(1).Height;
        info(1).BitDepth;
        CFP_img = zeros(info(1).Height,info(1).Width,length(info),'uint16');
        for i = 1:length(info)
            CFP_img(:,:,i) = imread(filename_CFP,i);
        end

        CFP_img_wcm = CFP_img;
        CFP_img_wcm(~wc_mask) = 0;

        savePath = sprintf('/Users/yijiachen/Documents/B_Cell_project/HVN037/images/HVN037_live_nucleus_mask/CFP_wholecell_masked/s%d/CFP_img_wc_masked_%d.tiff', s, t);
        imwrite(CFP_img_wcm, savePath);


        CFP_stats = regionprops("table",wc_mask,CFP_img_wcm,'Area',"Centroid",'MeanIntensity',"MajorAxisLength","MinorAxisLength",'PixelValues') ;
        CFP_stats = CFP_stats(CFP_stats.Area>=100 & CFP_stats.Area <= 800,:);
        diameters = mean([CFP_stats.MajorAxisLength CFP_stats.MinorAxisLength],2);
        radii_CFP = diameters/2;
        fn=['s', num2str(s)];

        for row_num = 1:size(CFP_stats,1)
            CFP_wcm_all_t.(fn).('wc_medianint')(row_num,time_point)=median(cell2mat(CFP_stats.PixelValues(row_num,1)));
            CFP_wcm_all_t.(fn).('wc_sumint')(row_num,time_point)=sum(cell2mat(CFP_stats.PixelValues(row_num,1)));
            CFP_wcm_all_t.(fn).('wc_meanint')(row_num,time_point)=CFP_stats.MeanIntensity(row_num,1);
            CFP_wcm_all_t.(fn).('wc_centroid')(row_num,((2*(time_point-1)+1):(2*(time_point-1)+1)+1))=CFP_stats.Centroid(row_num,1:2);
            CFP_wcm_all_t.(fn).('radius')(row_num,time_point)=radii_CFP(row_num,1);
            CFP_wcm_all_t.(fn).('area')(row_num,time_point)=CFP_stats.Area(row_num,1);
        end

        %% load YFP image
        filename_YFP = sprintf('HVN037a_NoStim_TFP_Venus_Cherry_IkBeKO_40XOil_RedBF-02_scene%d_t%02dphiYFP_2_ORG.tif', s, t);

        info = imfinfo(filename_YFP);
        length(info);
        info(1).Width;
        info(1).Height;
        info(1).BitDepth;
        YFP_img = zeros(info(1).Height,info(1).Width,length(info),'uint16');
        for i = 1:length(info)
            YFP_img(:,:,i) = imread(filename_YFP,i);
        end

        YFP_img_wcm = YFP_img;
        YFP_img_wcm(~wc_mask) = 0;

        % Save the image as a TIFF file without reducing the resolution
        savePath = sprintf('/Users/yijiachen/Documents/B_Cell_project/HVN037/images/HVN037_live_nucleus_mask/YFP_wholecell_masked/s%d/YFP_img_wc_masked_%d.tiff', s, t);
        imwrite(YFP_img_wcm, savePath);


        YFP_stats = regionprops("table",wc_mask,YFP_img_wcm,'Area',"Centroid",'MeanIntensity',"MajorAxisLength","MinorAxisLength",'PixelValues') ;
        YFP_stats = YFP_stats(YFP_stats.Area>=100 & YFP_stats.Area <= 800,:);
        diameters = mean([YFP_stats.MajorAxisLength YFP_stats.MinorAxisLength],2);
        radii_YFP = diameters/2;

        for row_num_2 = 1:size(YFP_stats,1)
            YFP_wcm_all_t.(fn).('wc_medianint')(row_num_2,time_point)=median(cell2mat(YFP_stats.PixelValues(row_num_2,1)));
            YFP_wcm_all_t.(fn).('wc_sumint')(row_num_2,time_point)=sum(cell2mat(YFP_stats.PixelValues(row_num_2,1)));
            YFP_wcm_all_t.(fn).('wc_meanint')(row_num_2,time_point)=YFP_stats.MeanIntensity(row_num_2,1);
            YFP_wcm_all_t.(fn).('wc_centroid')(row_num_2,((2*(time_point-1)+1):(2*(time_point-1)+1)+1))=YFP_stats.Centroid(row_num_2,1:2);
            YFP_wcm_all_t.(fn).('radius')(row_num_2,time_point)=radii_YFP(row_num_2,1);
            YFP_wcm_all_t.(fn).('area')(row_num,time_point)=YFP_stats.Area(row_num,1);
        end

        time_point=time_point+1;

% % %         %% visualization
% % %         figure
% % %         subplot(2,2,1);
% % %         imagesc(CFP_img_wcm)
% % %         axis off;
% % %         colormap gray
% % %         hold on
% % %         scatter(CFP_stats.Centroid(:, 1), CFP_stats.Centroid(:, 2), 'r.', 'SizeData', 20); % Adjust 'SizeData' as needed
% % %         % Plot the mean intensity values
% % %         for i = 1:size(CFP_stats,1)
% % %             x = CFP_stats.Centroid(i, 1); % X-coordinate of the centroid
% % %             y = CFP_stats.Centroid(i, 2); % Y-coordinate of the centroid
% % %             pixelvalue_CFP = median(cell2mat(CFP_stats.PixelValues(i))); % Mean intensity value
% % %             if pixelvalue_CFP >=4500
% % %                 % Add the mean intensity value as text at the centroid
% % %                 text(x, y+15, num2str(pixelvalue_CFP), 'Color', 'yellow', 'FontSize', 10, 'FontWeight', 'bold');
% % %                 pixelvalue_CFP
% % %             end
% % %         end
% % % 
% % %         subplot(2,2,3);
% % %         [CFP_img_display,map1]=imread(filename_CFP);
% % %         imshow(CFP_img_display,map1)
% % % 
% % %         subplot(2,2,2);
% % %         imagesc(YFP_img_wcm)
% % %         axis off;
% % %         colormap gray
% % %         hold on
% % %         scatter(YFP_stats.Centroid(:, 1), YFP_stats.Centroid(:, 2), 'r.', 'SizeData', 20); % Adjust 'SizeData' as needed
% % %         % Plot the mean intensity values
% % %         for i = 1:size(YFP_stats,1)
% % %             x = YFP_stats.Centroid(i, 1); % X-coordinate of the centroid
% % %             y = YFP_stats.Centroid(i, 2); % Y-coordinate of the centroid
% % %             pixelvalue_YFP = median(cell2mat(YFP_stats.PixelValues(i))); % Mean intensity value
% % %             % Add the mean intensity value as text at the centroid
% % %             if pixelvalue_YFP >=1000
% % %                 text(x, y+15, num2str(pixelvalue_YFP), 'Color', 'yellow', 'FontSize', 10, 'FontWeight', 'bold');
% % %             end
% % %         end
% % %         subplot(2,2,4);
% % %         [YFP_img_display,map2]=imread(filename_YFP);
% % %         imshow(YFP_img_display,map2)


    end

end


            save('/Users/yijiachen/Documents/B_Cell_project/HVN037/data_saved/HVN037_CFP_wholecell_mask_value.mat','CFP_wcm_all_t')
            save('/Users/yijiachen/Documents/B_Cell_project/HVN037/data_saved/HVN037_YFP_wholecell_mask_value.mat','YFP_wcm_all_t')
