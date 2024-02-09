for s = 8
    %1:9

    for t = 145

        %% APC death mask
        path = sprintf('/Users/yijiachen/Documents/B_Cell_project/panel_cd/HVN032_Fluo_9_scenes/HVN032_Fluo_ALL_s%d', s);
        addpath(path);
        filename_APC = ['HVN032_Fluo_ALL_t', num2str(t,'%03.f'), 'c2_ORG', '.tif'];
        info = imfinfo(filename_APC);
        length(info);
        info(1).Width;
        info(1).Height;
        info(1).BitDepth;
        APC_img = zeros(info(1).Height,info(1).Width,length(info),'uint16');
        for i = 1:length(info)
            APC_img(:,:,i) = imread(filename_APC,i);
        end

%                 % histogram of the image
%         % Compute the histogram
% histogram = imhist(APC_img);
% 
% % Plot the histogram
% figure;
% bar(histogram);
% title('APC_img Histogram');
% xlabel('Pixel Value');
% ylabel('Frequency');
% xlim([0 50])


        figure
        imagesc(APC_img) %default set the scale range to [min, max]
        truesize([500,500])
        colormap gray
        colorbar
        title("original image APC")

        % Otsu original
        APC_death_mask = imbinarize(uint16(APC_img),'global');
        figure
        imagesc(APC_death_mask)
        truesize([500 500]);
        colormap gray
        colorbar
        title("APC death mask otsu thresholding")



% %         % Adaptive original
%         APC_death_mask = imbinarize(uint16(APC_img),'adaptive');
%         figure
%         imagesc(APC_death_mask)
%         truesize([500 500]);
%         colormap gray
%         colorbar
%         title("APC death mask adaptive thresholding")
% % 
%         APC_death_mask= imopen(APC_death_mask, strel('disk', 3));
%         figure
%         imshow(APC_death_mask)
%         colormap gray
%         axis on
%         title("openned image APC")

        %watershed segmentation
        D = bwdist(~APC_death_mask);
        D = imcomplement(D);
        minima = imextendedmin(D,1);
        D = imimposemin(D,minima);
        W = watershed(D);

        APC_death_mask(W==0)=0;
        figure
        imagesc(APC_death_mask)
        truesize([500 500]);
        colormap gray
        colorbar
        title('APC death mask watershed')

        % figure
        % label = bwlabel(APC_death_mask);
        % s = regionprops(label, 'Centroid');
        % imshow(APC_death_mask)
        % hold on
        % for k = 1:numel(s)
        %     c = s(k).Centroid;
        %     text(c(1), c(2), sprintf('%d', k), ...
        %         'HorizontalAlignment', 'center', ...
        %         'VerticalAlignment', 'middle','Color','red');
        % end
        % hold off

        %% mCherry Nuclear Mask
        filename_mCherry = ['HVN032_Fluo_ALL_t', num2str(t,'%03.f'), 'c3_ORG', '.tif'];
        info = imfinfo(filename_mCherry);
        length(info);
        info(1).Width;
        info(1).Height;
        info(1).BitDepth;
        mCherry_img = zeros(info(1).Height,info(1).Width,length(info),'uint16');
        for i = 1:length(info)
            mCherry_img(:,:,i) = imread(filename_mCherry,i);
        end

        figure
        imagesc(mCherry_img, [400, 800]) %default set the scale range to [min, max]
        truesize([500,500])
        colormap gray
        colorbar
        title("original image mCherry")

        %Adaptive original image
        adaptive_it_ori = imbinarize(uint16(mCherry_img),'adaptive');
        figure
        imshow(adaptive_it_ori)
        title("adaptive thresholding original")

        % %erode
        % adaptive_it_ori_erode = imerode(adaptive_it_ori,strel('disk', 2));
        % figure
        % imshow(adaptive_it_ori_erode)
        % title("adaptive thresholding original erosion")
        %
        % adaptive_it_ori_dilate = imdilate(adaptive_it_ori_erode,strel('disk', 2));
        % figure
        % imshow(adaptive_it_ori_dilate)
        % title("adaptive thresholding original dilation")

        % %open adaptive
        % mCherry_nuclear_mask= imopen(adaptive_it_ori, strel('disk', 3));
        % figure
        % imshow(mCherry_nuclear_mask)
        % imagesc(mCherry_nuclear_mask)
        % truesize([500 500]);
        % colormap gray
        % colorbar
        % title("mCherry nuclear mask")

        mCherry_nuclear_mask= imopen(adaptive_it_ori, strel('disk', 3));
        figure
        imshow(mCherry_nuclear_mask)
        colormap gray
        axis on
        title("openned image mCherry")
        %watershed segmentation
        D = bwdist(~mCherry_nuclear_mask);
        D = imcomplement(D);
        minima = imextendedmin(D,1);
        D = imimposemin(D,minima);
        W = watershed(D);

        % figure
        mCherry_nuclear_mask(W==0)=0;
        imagesc(mCherry_nuclear_mask)
        truesize([500 500]);
        colormap gray
        colorbar
        title('mcherry nuclear mask watershed')

        % figure
        % label = bwlabel(mCherry_nuclear_mask);
        % s = regionprops(label, 'Centroid');
        % imshow(mCherry_nuclear_mask)
        % hold on
        % for k = 1:numel(s)
        %     c = s(k).Centroid;
        %     text(c(1), c(2), sprintf('%d', k), ...
        %         'HorizontalAlignment', 'center', ...
        %         'VerticalAlignment', 'middle','Color','red');
        % end
        % hold off

        %% Live Mask = xor (mCherry nuclear mask, APC death mask)
        live_mask_xor = xor(mCherry_nuclear_mask,APC_death_mask);
        % figure
        % imagesc(live_mask_xor)
        % truesize([500 500]);
        % colormap gray
        % colorbar
        % title('live mask xor')

        %eliminate the ring structure
        no_holes = bwpropfilt(live_mask_xor,'EulerNumber',[1 1]);
        % figure
        % imagesc(no_holes)
        % truesize([500 500]);
        % colormap gray
        % colorbar
        % title('live mask xor no holes')

        %further eliminate the ring structure

        final_live_mask = imopen (no_holes,strel('disk', 3));

        figure
        imagesc(final_live_mask)
        truesize([500,500])
        colormap gray
        colorbar
        title("final live mask")

        %% label the cells
        figure
        label = bwlabel(final_live_mask);
        s_label = regionprops(label, 'Centroid');
        imshow(final_live_mask)
        hold on
        for k = 1:numel(s_label)
            c = s_label(k).Centroid;
            text(c(1), c(2), sprintf('%d', k), ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle','Color','red');
        end
        hold off

        %% load CFP image
        filename_CFP = ['HVN032_Fluo_ALL_t', num2str(t,'%03.f'), 'c5_ORG', '.tif'];
        info = imfinfo(filename_CFP);
        length(info);
        length(info);
        info(1).Width;
        info(1).Height;
        info(1).BitDepth;
        CFP_img = zeros(info(1).Height,info(1).Width,length(info),'uint16');
        for i = 1:length(info)
            CFP_img(:,:,i) = imread(filename_CFP,i);
        end
        % figure
        % imagesc(CFP_img) %default set the scale range to [min, max]
        % truesize([500,500])
        % colormap gray
        % colorbar
        % title("original image CFP")

        CFP_img_mask = CFP_img;
        CFP_img_mask(~final_live_mask) = 0;
        % CFP_img_mask = CFP_img.*final_live_mask;

        % figure
        % imagesc(CFP_img_mask) %default set the scale range to [min, max]
        % truesize([500,500])
        % colormap gray
        % colorbar
        % title("CFP image live mask")

        label_CFP = bwlabel(CFP_img_mask);
        s_CFP = regionprops(label_CFP, 'Centroid');
        CFP_value = regionprops("table",label_CFP,CFP_img_mask,'MeanIntensity') ;
        figure
        imagesc(CFP_img_mask)
        truesize([500,500])
        colormap gray
        colorbar
        hold on
        for k = 1:numel(s_CFP)
            c = s_CFP(k).Centroid;
            %     mean = CFP_value.MeanIntensity(k);
            text(c(1), c(2), sprintf('%d', k), ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle','Color','red');
        end
        hold off
        title("CFP image live mask labeled")

        %% load YFP image
        filename_YFP = ['HVN032_Fluo_ALL_t', num2str(t,'%03.f'), 'c4_ORG', '.tif'];

        info = imfinfo(filename_YFP);length(info);
        info(1).Width;
        info(1).Height;
        info(1).BitDepth;
        YFP_img = zeros(info(1).Height,info(1).Width,length(info),'uint16');
        for i = 1:length(info)
            YFP_img(:,:,i) = imread(filename_YFP,i);
        end
        % figure
        % imagesc(YFP_img) %default set the scale range to [min, max]
        % truesize([500,500])
        % colormap gray
        % colorbar
        % title("original image YFP")

        YFP_img_mask = YFP_img;
        YFP_img_mask(~final_live_mask) = 0;
        % YFP_img_mask = YFP_img.*final_live_mask;
        % figure
        % imagesc(YFP_img_mask) %default set the scale range to [min, max]
        % truesize([500,500])
        % colormap gray
        % colorbar
        % title("YFP image live mask")

        label_YFP = bwlabel(YFP_img_mask);
        s_YFP = regionprops(label_YFP, 'Centroid');
        YFP_value = regionprops("table",label_YFP,YFP_img_mask,'MeanIntensity') ;
        figure
        imagesc(YFP_img_mask)
        truesize([500,500])
        colormap gray
        colorbar
        hold on
        for k = 1:numel(s_YFP)
            c = s_YFP(k).Centroid;
            %     mean = YFP_value.MeanIntensity(k);
            text(c(1), c(2), sprintf('%d', k), ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle','Color','red');
        end
        hold off
        title("YFP image live mask labeled")
    end
end