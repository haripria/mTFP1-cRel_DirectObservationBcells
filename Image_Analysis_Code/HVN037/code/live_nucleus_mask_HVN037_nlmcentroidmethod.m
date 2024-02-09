for s = [1,2,5,6,7,8]
    %1:9

    for t = 1:92
        %1:30:150
        %1:150

        path = sprintf('/Volumes/Data/Mark/HVN037a_NoStim_TFP_Venus_Cherry_IkBeKO_40XOil_RedBF-02_scene%d', s);
        addpath(path);
%% No APC for HVN037
%         filename_APC = ['HVN037_Fluo_ALL_t', num2str(t,'%03.f'), 'c2_ORG', '.tif'];
%         info = imfinfo(filename_APC);
%         length(info);
%         info(1).Width;
%         info(1).Height;
%         info(1).BitDepth;
%         APC_img = zeros(info(1).Height,info(1).Width,length(info),'uint16');
%         for i = 1:length(info)
%             APC_img(:,:,i) = imread(filename_APC,i);
%         end
% 
% 
%         %         % Otsu original
%         %         APC_death_mask = imbinarize(uint16(APC_img),'global');
%         % Adaptive original
%         APC_death_mask = imbinarize(uint16(APC_img),'adaptive');
%         APC_death_mask= imopen(APC_death_mask, strel('disk', 3));
% 
%         savePath = sprintf('/Users/yijiachen/Documents/B_Cell_project/HVN037/images/HVN037_live_nucleus_mask/APC_death_mask/s%d/APC_death_mask_%d.tiff',s, t);
%         imwrite(APC_death_mask, savePath);

        %% mCherry Nuclear Mask        
%         filename_mCherry = ['HVN037a_NoStim_TFP_Venus_Cherry_IkBeKO_40XOil_RedBF-02_scene1_t', num2str(t,'%02.f'), 'mCherry_1_ORG', '.tif'];
        filename_mCherry = sprintf('HVN037a_NoStim_TFP_Venus_Cherry_IkBeKO_40XOil_RedBF-02_scene%d_t%02dmCherry_1_ORG.tif', s, t);

        info = imfinfo(filename_mCherry);
        length(info);
        info(1).Width;
        info(1).Height;
        info(1).BitDepth;
        mCherry_img = zeros(info(1).Height,info(1).Width,length(info),'uint16');
        for i = 1:length(info)
            mCherry_img(:,:,i) = imread(filename_mCherry,i);
        end

        %Adaptive original image
        adaptive_it_ori = imbinarize(uint16(mCherry_img),'adaptive');

        adaptive_it_ori_erode = imerode(adaptive_it_ori,strel('disk', 2));
        %open adaptive
        mCherry_nuclear_mask= imopen(adaptive_it_ori_erode, strel('disk', 1));

        %         mCherry_nuclear_mask= imopen(adaptive_it_ori, strel('disk', 3));

        %watershed segmentation
        D = bwdist(~mCherry_nuclear_mask);
        D = imcomplement(D);
        minima = imextendedmin(D,1);
        D = imimposemin(D,minima);
        W = watershed(D);

        mCherry_nuclear_mask(W==0)=0;


        savePath = sprintf('/Users/yijiachen/Documents/B_Cell_project/HVN037/images/HVN037_live_nucleus_mask/mCherry_nucleus_mask/s%d/mCherry_nucleus_mask_%d.tiff',s, t);
        imwrite(mCherry_nuclear_mask, savePath);


        %% Live Mask (NA because no APC)
%         label_mCherry = bwlabel(mCherry_nuclear_mask);
%         mCherry_value = regionprops("table",label_mCherry,mCherry_nuclear_mask,'Centroid') ;

% 
%         % Get a list of unique labels
%         uniqueLabels = unique(label_mCherry);
% 
%         % Find the index of label 0 in uniqueLabels
%         zeroLabelIndex = find(uniqueLabels == 0);
% 
%         % Exclude label 0 if it exists in uniqueLabels
%         if ~isempty(zeroLabelIndex)
%             uniqueLabels(zeroLabelIndex) = [];
%         end
% 
%         % Initialize a binary mask with zeros
%         final_live_mask = zeros(size(label_mCherry), 'logical');
%         mCherry_value_array = table2array(mCherry_value);
%         % Create a combined binary mask
%         for i = 1:numel(uniqueLabels)
%             if APC_death_mask(round(mCherry_value_array(i,2)),round(mCherry_value_array(i,1)))==0
%                 label = uniqueLabels(i);
%                 binaryMask = label_mCherry == label;
%                 final_live_mask = final_live_mask | binaryMask;
%             end
%         end
% 
% 
%         savePath = sprintf('/Users/yijiachen/Documents/B_Cell_project/HVN037/images/HVN037_live_nucleus_mask/live_nucleus_mask/s%d/live_nucleus_mask_%d.tiff',s, t);
%         imwrite(final_live_mask, savePath);

% % %         %% visualization
% % %                         figure
% % %                         imagesc(mCherry_nuclear_mask)
% % %                         truesize([500 500]);
% % %                         colormap gray
% % %                         colorbar
% % %                         title("mCherry adaptive thresholding")
% % %                         hold on
% % %         
% % %         
% % %                         % Plot centroids as red dots on the image
% % %                         centroid_y = mCherry_value.Centroid(:, 1);
% % %                         centroid_x = mCherry_value.Centroid(:, 2);
% % %                         scatter(centroid_y, centroid_x, 'r.', 'SizeData', 100); % Adjust 'SizeData' as needed
% % %         
% % %         
% % % % %                         figure
% % % % %                         imagesc(APC_death_mask)
% % % % %                         truesize([500 500]);
% % % % %                         colormap gray
% % % % %                         colorbar
% % % % %                         title("APC death mask adaptive thresholding")
% % % % %                         hold on
% % % % %                         scatter(centroid_y, centroid_x, 'r.', 'SizeData', 100); % Adjust 'SizeData' as needed
% % % % %         
% % % % %                         figure
% % % % %                         imagesc(final_live_mask)
% % % % %                         truesize([500 500]);
% % % % %                         colormap gray
% % % % %                         colorbar
% % % % %                         title("final_live_mask")
% % % % %                         hold on
% % % % %                         scatter(centroid_y, centroid_x, 'r.',
% % % % %         'SizeData', 100); % Adjust 'SizeData' as needed
% % %                 % Raw Images
% % %              
% % % %                 filename_BF = ['HVN037a_NoStim_TFP_Venus_Cherry_IkBeKO_40XOil_RedBF-02_scene1_t', num2str(t,'%03.f'), '_ORG', '.tif'];
% % %                 filename_BF = sprintf('HVN037a_NoStim_TFP_Venus_Cherry_IkBeKO_40XOil_RedBF-02_scene%d_t%03d_ORG.tif', s, t);
% % %                 [BF_img_display,map3]=imread(filename_BF);
% % % % %                 [APC_img_display,map1]=imread(filename_APC);
% % %                 [mCherry_img_display,map2]=imread(filename_mCherry);
% % %                 figure
% % % % %                 subplot(1,3,2), imshow(APC_img_display,map1)
% % %                 subplot(1,3,1), imshow(mCherry_img_display,map2)
% % %                 subplot(1,3,3), imshow(BF_img_display,map3)

    end
end
