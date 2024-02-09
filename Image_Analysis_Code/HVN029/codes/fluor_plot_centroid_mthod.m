

for cell_num = 17
    %:17


    load('/Users/yijiachen/Documents/B_Cell_project/HVN029/trackmate_data/HVN029_BF_ALL_s2_t0001_t3000_ORG_Concat_trial2_trackdata.mat'); % Load the centroid data

    CFP_MeanIntensity = [];

    frame_tot = alltrackData{cell_num}.Frame;


    centroid_tot = [alltrackData{cell_num}.X alltrackData{cell_num}.Y];

    for r = 1:size(frame_tot,1)
        frame_num = frame_tot(r,1)-1;
        if mod(frame_num - 1, 20) == 0
            t = (frame_num - 1) / 20 + 1;


            addpath '/Users/yijiachen/Documents/B_Cell_project/HVN029/images/HVN029_live_nucleus_mask/CFP_masked/s2'

            filename_CFP = sprintf('CFP_img_masked_%d.tiff', t);
            info = imfinfo(filename_CFP);
            length(info);
            info(1).Width;
            info(1).Height;
            info(1).BitDepth;
            CFP_img = zeros(info(1).Height,info(1).Width,length(info),'uint16');
            for read = 1:length(info)
                CFP_img(:,:,read) = imread(filename_CFP,read);
            end
            label_CFP = bwlabel(CFP_img);

            CFP_stats = regionprops("table",label_CFP,CFP_img,"Centroid","MajorAxisLength","MinorAxisLength",'MeanIntensity') ;
            centroid_nucl_1 = CFP_stats.Centroid;
            diameters_1 = mean([CFP_stats.MajorAxisLength CFP_stats.MinorAxisLength],2);
            radii_1 = diameters_1/2;
            meanintensity_1 = CFP_stats.MeanIntensity;
            centroid_diff_1 = zeros(size(centroid_nucl_1,1),1);

            for c1 = 1:size(centroid_nucl_1,1)
                centroid_diff_1(c1,1) = ((centroid_tot(r,1)-centroid_nucl_1(c1,1))^2 + (centroid_tot(r,2)-centroid_nucl_1(c1,2))^2)^(1/2);
            end


            CFP_MeanIntensity(t,1)=frame_num;
            if ~isempty(find(centroid_diff_1(:,1)<= radii_1(:,1), 1))
                CFP_MeanIntensity(t,2)=meanintensity_1(centroid_diff_1(:,1)<= radii_1(:,1),1);

            else
                CFP_MeanIntensity(t,2)=NaN;

            end

            addpath '/Users/yijiachen/Documents/B_Cell_project/HVN029/images/HVN029_live_nucleus_mask/YFP_masked/s2'
            filename_YFP = sprintf('YFP_img_masked_%d.tiff', t);
            info = imfinfo(filename_YFP);
            length(info);
            info(1).Width;
            info(1).Height;
            info(1).BitDepth;
            YFP_img = zeros(info(1).Height,info(1).Width,length(info),'uint16');
            for read = 1:length(info)
                YFP_img(:,:,read) = imread(filename_YFP,read);
            end
            label_YFP = bwlabel(YFP_img);

            YFP_stats = regionprops("table",label_YFP,YFP_img,"Centroid","MajorAxisLength","MinorAxisLength",'MeanIntensity') ;
            centroid_nucl_2 = YFP_stats.Centroid;
            diameters_2 = mean([YFP_stats.MajorAxisLength YFP_stats.MinorAxisLength],2);
            radii_2 = diameters_2/2;
            meanintensity_2 = YFP_stats.MeanIntensity;
            centroid_diff_2 = zeros(size(centroid_nucl_2,1),1);

            for c2 = 1:size(centroid_nucl_2,1)
                centroid_diff_2(c2,1) = ((centroid_tot(r,1)-centroid_nucl_2(c2,1))^2 + (centroid_tot(r,2)-centroid_nucl_2(c2,2))^2)^(1/2);
            end


            YFP_MeanIntensity(t,1)=frame_num;
            if ~isempty(find(centroid_diff_2(:,1)<= radii_2(:,1), 1))
                YFP_MeanIntensity(t,2)=meanintensity_2(centroid_diff_2(:,1)<= radii_2(:,1),1);

            else
                YFP_MeanIntensity(t,2)=NaN;

            end
% % %             %% visualization
% % %             centroid_benchmark = figure;
% % %             subplot(1,2,1)
% % %             imshow(CFP_img, []);
% % %             hold on;
% % %             scatter(centroid_tot(t, 1), centroid_tot(t, 2),1,".");
% % %             hold off;
% % %             title(sprintf('t = %d CFP Image',t));
% % %             subplot(1,2,2)
% % %             imshow(YFP_img, []);
% % %             hold on;
% % %             scatter(centroid_tot(t, 1), centroid_tot(t, 2),1,".");
% % %             hold off;
% % %             title('YFP Image');
% % %             saveas(centroid_benchmark,sprintf('/Users/yijiachen/Documents/B_Cell_project/HVN029/images/benchmark/s2/c%d/centroid_benchmark_c%d_t%d.png',cell_num,cell_num,t))
% % %             close all
        end
    end
    YFP_CFP_nuclear_plot = figure;

    valid_indices_CFP = ~isnan(CFP_MeanIntensity(:, 2)) & CFP_MeanIntensity(:, 2) ~= 0;

    plot(CFP_MeanIntensity(valid_indices_CFP, 1), CFP_MeanIntensity(valid_indices_CFP, 2), 'LineWidth', 2);

    hold on

    valid_indices_YFP = ~isnan(YFP_MeanIntensity(:, 2)) & YFP_MeanIntensity(:, 2) ~= 0;

    plot(YFP_MeanIntensity(valid_indices_YFP, 1), 3.5 * YFP_MeanIntensity(valid_indices_YFP, 2), 'LineWidth', 2);


    ylim([2300 3200])
    xlim([0 1000])
    xticks([1 120 240 360 480 600 720 840 960])
    xticklabels({'0','2','4','6','8','10','12','14','16'})
    title('YFP/CFP nuclear mean intensity')
    xlabel('time (h)')
    ylabel('YFP/CFP nuclear sum intensity')
    legend ('CFP','YFP')
    saveas(YFP_CFP_nuclear_plot,sprintf('/Users/yijiachen/Documents/B_Cell_project/HVN029/images/s2/fluor_plot_centroid_method_c%d.png',cell_num))

    clear

end
