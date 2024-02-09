for s = 1:9
    time_point=1;
    for t = [1 2 3 4 7 13 19 25 49 97 145]
        path = sprintf('/Users/yijiachen/Documents/B_Cell_project/panel_cd/HVN032_live_nucleus_mask/APC_death_mask/s%d',s);
        addpath(path);
        filename_APC_death_mask = ['APC_death_mask_', num2str(t,'%01.f'), '.tiff'];
        info = imfinfo(filename_APC_death_mask);
        length(info);
        info(1).Width;
        info(1).Height;
        info(1).BitDepth;
        APC_death_mask = zeros(info(1).Height,info(1).Width,length(info),'uint16');
        for i = 1:length(info)
            APC_death_mask(:,:,i) = imread(filename_APC_death_mask,i);
        end

        path = sprintf('/Users/yijiachen/Documents/B_Cell_project/panel_cd/HVN032_BF_mask/s%d', s);
        addpath(path);
        filename_cellpose_mask = ['HVN032_Fluo_ALL_t', num2str(t,'%03.f'), '_cp_masks', '.png'];

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
            if APC_death_mask(round(cellpose_value_array(i,2)),round(cellpose_value_array(i,1)))==0
                label = uniqueLabels(i);
                nonbinaryMask = uint16(cellpose_mask == label)*label;
                cellpose_mask_final = cellpose_mask_final + nonbinaryMask;
            end
        end

        savePath = sprintf('/Users/yijiachen/Documents/B_Cell_project/panel_cd/HVN032_live_nucleus_mask/cellpose_wholecell_masked/s%d/cellpose_wholecell_mask_%d.tiff',s, t);
        imwrite(cellpose_mask_final, savePath);

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

        CFP_img_mask = CFP_img;
        CFP_img_mask(~cellpose_mask_final) = 0;

        % Save the image as a TIFF file without reducing the resolution
        savePath = sprintf('/Users/yijiachen/Documents/B_Cell_project/panel_cd/HVN032_live_nucleus_mask/CFP_wholecell_masked/s%d/CFP_wholecell_masked_%d.tiff', s, t);
        imwrite(CFP_img_mask, savePath);

        CFP_stats = regionprops("table",cellpose_mask_final,CFP_img_mask,'Area','Centroid','PixelValues') ;



        %% intensity data recorded
        field_name = ['s', num2str(s)];
        fn=['s', num2str(s)];
        for row_num = 1:size(CFP_stats,1)
            CFP_wholecell_all_t.(fn).('wc_medianint')(row_num,time_point)=median(cell2mat(CFP_stats.PixelValues(row_num,1)));
            CFP_wholecell_all_t.(fn).('wc_sumint')(row_num,time_point)=sum(cell2mat(CFP_stats.PixelValues(row_num,1)));
            CFP_wholecell_all_t.(fn).('wc_area')(row_num,time_point)=CFP_stats.Area(row_num,1);
            CFP_wholecell_all_t.(fn).('wc_centroid')(row_num,((2*(time_point-1)+1):(2*(time_point-1)+1)+1))=CFP_stats.Centroid(row_num,1:2);
        end
        time_point=time_point+1;

        % % %         %% visualization
        % % %         figure
        % % %         imagesc(cellpose_mask)
        % % %         truesize([500 500]);
        % % %         colormap gray
        % % %         colorbar
        % % %         title("mCherry adaptive thresholding")
        % % %         hold on
        % % %         scatter(centroid_y, centroid_x, 'r.', 'SizeData', 100); % Adjust 'SizeData' as needed
        % % %
        % % %         figure
        % % %         imagesc(APC_death_mask)
        % % %         truesize([500 500]);
        % % %         colormap gray
        % % %         colorbar
        % % %         title("APC death mask adaptive thresholding")
        % % %         hold on
        % % %         scatter(centroid_y, centroid_x, 'r.', 'SizeData', 100); % Adjust 'SizeData' as needed
        % % %
        % % %         %% visualization
        % % %         figure
        % % %         imagesc(cellpose_mask_final)
        % % %         truesize([500 500]);
        % % %         colormap gray
        % % %         colorbar
        % % %         title("final_cellpose_Mask")
        % % %         hold on
        % % %         scatter(centroid_y, centroid_x, 'r.', 'SizeData', 100); % Adjust 'SizeData' as needed
        % % %
        % % %
        % % %         %% visualization
        % % %         figure
        % % %         imagesc(CFP_img_mask)
        % % %         truesize([500 500]);
        % % %         colormap gray
        % % %         colorbar
        % % %         title("CFP_img_mask")
        % % %         hold on
        % % %         scatter(centroid_y, centroid_x, 'r.', 'SizeData', 100); % Adjust 'SizeData' as needed
        % % %         % Plot the mean intensity values
        % % %         for i = 1:size(CFP_stats,1)
        % % %             x = CFP_stats.Centroid(i, 1); % X-coordinate of the centroid
        % % %             y = CFP_stats.Centroid(i, 2); % Y-coordinate of the centroid
        % % % %             meanIntensity = floor(CFP_stats.MeanIntensity(i)); % Mean intensity value
        % % % %             % Add the mean intensity value as text at the centroid
        % % % %             text(x, y+15, num2str(meanIntensity), 'Color', 'yellow', 'FontSize', 10, 'FontWeight', 'bold');
        % % % %             areavalue = floor(CFP_stats.Area(i)); % Mean intensity value
        % % % %             % Add the mean intensity value as text at the centroid
        % % % %             text(x, y+15, num2str(areavalue), 'Color', 'yellow', 'FontSize', 10, 'FontWeight', 'bold');
        % % %             pixelvalue = sum(cell2mat(CFP_stats.PixelValues(i))); % Mean intensity value
        % % %             % Add the mean intensity value as text at the centroid
        % % %             text(x, y+15, num2str(pixelvalue), 'Color', 'yellow', 'FontSize', 10, 'FontWeight', 'bold');
        % % %         end


    end

        save('/Users/yijiachen/Documents/B_Cell_project/panel_cd/data_saved/HVN032_CFP_wholecell_value.mat','CFP_wholecell_all_t')
end