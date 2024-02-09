YFP_img_bg_mask_nonzero_mean=zeros(8,11);
YFP_img_bg_mask_nonzero_median=zeros(8,11);
for s = 1:8
    time_point=1;

    for t = [1 2 3 4 7 13 19 25 49 97 145]

        path = sprintf('/Users/yijiachen/Documents/B_Cell_project/panel_cd/HVN029_BF_mask/s%d', s);
        addpath(path);
        filename_cellpose_mask = ['HVN029_Fluo_ALL_t', num2str(t,'%03.f'), '_cp_masks', '.png'];

        cellpose_mask=imread(filename_cellpose_mask);

        cellpose_value = regionprops("table",cellpose_mask,cellpose_mask,'Centroid') ;
        % Plot centroids as red dots on the image
        centroid_y = cellpose_value.Centroid(:, 1);
        centroid_x = cellpose_value.Centroid(:, 2);

        % Get a list of unique labels
        uniqueLabels = unique(cellpose_mask);

        % Find the index of label 0 in uniqueLabels
        zeroLabelIndex = find(uniqueLabels == 0);

        % Exclude label 0 if it exists in uniqueLabels
        if ~isempty(zeroLabelIndex)
            uniqueLabels(zeroLabelIndex) = [];
        end

        % Initialize a binary mask with zeros
        cellpose_mask_final = zeros(size(cellpose_mask), 'uint16');
        cellpose_value_array = table2array(cellpose_value);
        % Create a combined binary mask
        for i = 1:numel(uniqueLabels)
            %             if APC_death_mask(round(cellpose_value_array(i,2)),round(cellpose_value_array(i,1)))==0
            label = uniqueLabels(i);
            nonbinaryMask = uint16(cellpose_mask == label)*label;
            cellpose_mask_final = cellpose_mask_final + nonbinaryMask;
            %             end
        end

        %         savePath = sprintf('/Users/yijiachen/Documents/B_Cell_project/panel_cd/HVN029_live_nucleus_mask/cellpose_wholecell_masked/s%d/cellpose_wholecell_mask_%d.tiff',s, t);
        %         imwrite(cellpose_mask_final, savePath);


        %% load YFP image
        path_2 = sprintf('/Users/yijiachen/Documents/B_Cell_project/panel_cd/HVN029_Fluo_8_scenes/HVN029_Fluo_ALL_s%d', s);
        addpath(path_2);

        filename_YFP = ['HVN029_Fluo_ALL_t', num2str(t,'%03.f'), 'c3_ORG', '.tif'];
        info = imfinfo(filename_YFP);
        length(info);
        info(1).Width;
        info(1).Height;
        info(1).BitDepth;
        YFP_img = zeros(info(1).Height,info(1).Width,length(info),'uint16');
        for i = 1:length(info)
            YFP_img(:,:,i) = imread(filename_YFP,i);
        end

        YFP_img_bg_mask = YFP_img;
        YFP_img_bg_mask(cellpose_mask_final~=0) = 0;


        % Assume the bright circle is centered and has a radius of half the minimum dimension of the image
        radius = 360;
        %min(size(image_gray)) / 2;
        center_x = 510;
        %size(image_gray, 2) / 2;
        center_y = 510;
        %size(image_gray, 1) / 2;

        % Create a grid of coordinates corresponding to the image dimensions
        [X, Y] = meshgrid(1:size(YFP_img_bg_mask, 2), 1:size(YFP_img_bg_mask, 1));

        % Calculate the distance of each point from the center
        distances = sqrt((X - center_x).^2 + (Y - center_y).^2);

        % Create a logical mask where points within the circle have value 1, others 0
        well_mask = distances <= radius;

        % Apply the mask to the image
        % Pixels outside the circle will be set to zero
        YFP_img_bg_mask(~well_mask) = 0;

        % Display the modified image
%         imshow(YFP_img_bg_mask);

        YFP_img_bg_mask_nonzero = YFP_img_bg_mask(YFP_img_bg_mask~=0);
        YFP_img_bg_mask_nonzero_mean(s,time_point) = mean(YFP_img_bg_mask_nonzero);
        YFP_img_bg_mask_nonzero_median(s,time_point) = median(YFP_img_bg_mask_nonzero);
        % Save the image as a TIFF file without reducing the resolution
        savePath = sprintf('/Users/yijiachen/Documents/B_Cell_project/panel_cd/HVN029_cellpose_bg_mask/s%d/YFP_cellpose_bg_mask_%d.tiff', s, t);
        imwrite(YFP_img_bg_mask, savePath);

   
        time_point = time_point+1;

    end

    %     save('/Users/yijiachen/Documents/B_Cell_project/panel_cd/data_saved/HVN029_YFP_wholecell_value.mat','YFP_wholecell_all_t')
end

bg_mean_HVN029 = mean(YFP_img_bg_mask_nonzero_mean,'all');
bg_median_HVN029 = median(YFP_img_bg_mask_nonzero_median,'all');
save('/Users/yijiachen/Documents/B_Cell_project/panel_cd/data_saved/HVN029_YFP_bg_mean_median.mat','bg_mean_HVN029','bg_median_HVN029')

