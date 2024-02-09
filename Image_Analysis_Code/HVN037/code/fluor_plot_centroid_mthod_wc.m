tot_cell_num = 1;
for s =[1,2,5,6,7,8]

    trackmate_filmname = sprintf('/Users/yijiachen/Documents/B_Cell_project/HVN037/trackmate_data/HVN037_BF_All_s%d_ORG_concat_trackdata.mat', s);
    load(trackmate_filmname)

    for cell_num =  1:size(alltrackData,2)
        CFP_MeanIntensity = [];
        CFP_MedianIntensity = [];
        CFP_SumIntensity = [];
        CFP_ConcIntensity = [];

        YFP_MeanIntensity = [];
        YFP_MedianIntensity = [];
        YFP_SumIntensity = [];
        YFP_ConcIntensity = [];


        frame_tot = alltrackData{cell_num}.Frame;

        centroid_tot = [alltrackData{cell_num}.X alltrackData{cell_num}.Y];

        % figure;

        for r = 1:size(frame_tot,1)
            frame_num = frame_tot(r,1) - 1;
            if mod(frame_num - 1, 10) == 0
                t = (frame_num - 1) / 10 + 1;
                path_CFP = sprintf('/Users/yijiachen/Documents/B_Cell_project/HVN037/images/HVN037_live_nucleus_mask/CFP_wholecell_masked/s%d', s);
                addpath(path_CFP);
                %             addpath '/Users/yijiachen/Documents/B_Cell_project/HVN037/images/HVN037_live_nucleus_mask/CFPwholecell_masked/s1'

                filename_CFP = sprintf('CFP_img_wc_masked_%d.tiff', t);
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

                CFP_stats = regionprops("table",label_CFP,CFP_img,'Area',"Centroid","MajorAxisLength","MinorAxisLength",'MeanIntensity','PixelValues') ;
                centroid_wc_1 = CFP_stats.Centroid;
                diameters_1 = mean([CFP_stats.MajorAxisLength CFP_stats.MinorAxisLength],2);
                radii_1 = diameters_1/2;
                meanintensity_1 = CFP_stats.MeanIntensity;
                area_1 = CFP_stats.Area;

                for row_num = 1:size(CFP_stats,1)
                    medianintensity_1(row_num,1) = median(cell2mat(CFP_stats.PixelValues(row_num,1)));
                    sumintensity_1(row_num,1) = sum(cell2mat(CFP_stats.PixelValues(row_num,1)));
                end

                concentration_1 = sumintensity_1 ./ area_1;

                centroid_diff_1 = zeros(size(centroid_wc_1,1),1);


                for c1 = 1:size(centroid_wc_1,1)
                    centroid_diff_1(c1,1) = ((centroid_tot(r,1)-centroid_wc_1(c1,1))^2 + (centroid_tot(r,2)-centroid_wc_1(c1,2))^2)^(1/2);
                end


                CFP_MeanIntensity(t,1)=frame_num;
                CFP_MedianIntensity(t,1)=frame_num;
                CFP_SumIntensity(t,1)=frame_num;
                CFP_ConcIntensity(t,1)=frame_num;
                if ~isempty(find(centroid_diff_1(:,1)<= radii_1(:,1), 1))
                    %                     CFP_MeanIntensity(t,2)=meanintensity_1(centroid_diff_1(:,1)<= radii_1(:,1),1);
                    %                     CFP_MedianIntensity(t,2)=medianintensity_1(centroid_diff_1(:,1)<= radii_1(:,1),1);
                    %                     CFP_SumIntensity(t,2)=sumintensity_1(centroid_diff_1(:,1)<= radii_1(:,1),1);
                    %                     CFP_ConcIntensity(t,2)=concentration_1(centroid_diff_1(:,1)<= radii_1(:,1),1);

                    CFP_MeanIntensity(t,2)=meanintensity_1(centroid_diff_1(:,1)== min(centroid_diff_1(:,1)),1);
                    CFP_MedianIntensity(t,2)=medianintensity_1(centroid_diff_1(:,1)== min(centroid_diff_1(:,1)),1);
                    CFP_SumIntensity(t,2)=sumintensity_1(centroid_diff_1(:,1)== min(centroid_diff_1(:,1)),1);
                    CFP_ConcIntensity(t,2)=concentration_1(centroid_diff_1(:,1)== min(centroid_diff_1(:,1)),1);
                else
                    CFP_MeanIntensity(t,2)=NaN;
                    CFP_MedianIntensity(t,2)=NaN;
                    CFP_SumIntensity(t,2)=NaN;
                    CFP_ConcIntensity(t,2)=NaN;

                end

                %             addpath '/Users/yijiachen/Documents/B_Cell_project/HVN037/images/HVN037_live_nucleus_mask/YFP_masked/s1'
                path_YFP = sprintf('/Users/yijiachen/Documents/B_Cell_project/HVN037/images/HVN037_live_nucleus_mask/YFP_wholecell_masked/s%d', s);
                addpath(path_YFP);
                filename_YFP = sprintf('YFP_img_wc_masked_%d.tiff', t);
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
                YFP_stats = regionprops("table",label_YFP,YFP_img,'Area',"Centroid","MajorAxisLength","MinorAxisLength",'MeanIntensity','PixelValues') ;

                %                 YFP_stats = regionprops("table",label_YFP,YFP_img,"Centroid","MajorAxisLength","MinorAxisLength",'MeanIntensity') ;
                centroid_wc_2 = YFP_stats.Centroid;
                diameters_2 = mean([YFP_stats.MajorAxisLength YFP_stats.MinorAxisLength],2);
                radii_2 = diameters_2/2;
                meanintensity_2 = YFP_stats.MeanIntensity;
                area_2 = YFP_stats.Area;
                for row_num = 1:size(YFP_stats,1)
                    medianintensity_2(row_num,1) = median(cell2mat(YFP_stats.PixelValues(row_num,1)));
                    sumintensity_2(row_num,1) = sum(cell2mat(YFP_stats.PixelValues(row_num,1)));
                end

                concentration_2 = sumintensity_2 ./ area_2;
                centroid_diff_2 = zeros(size(centroid_wc_2,1),1);

                for c2 = 1:size(centroid_wc_2,1)
                    centroid_diff_2(c2,1) = ((centroid_tot(r,1)-centroid_wc_2(c2,1))^2 + (centroid_tot(r,2)-centroid_wc_2(c2,2))^2)^(1/2);
                end


                YFP_MeanIntensity(t,1)=frame_num;
                YFP_MedianIntensity(t,1)=frame_num;
                YFP_SumIntensity(t,1)=frame_num;
                YFP_ConcIntensity(t,1)=frame_num;
                if ~isempty(find(centroid_diff_2(:,1)<=  radii_2(:,1), 1))
                    %                     YFP_MeanIntensity(t,2)=meanintensity_2(centroid_diff_2(:,1)<= radii_2(:,1),1);
                    %                     YFP_MedianIntensity(t,2)=medianintensity_2(centroid_diff_2(:,1)<= radii_2(:,1),1);
                    %                     YFP_SumIntensity(t,2)=sumintensity_2(centroid_diff_2(:,1)<= radii_2(:,1),1);
                    %                     YFP_ConcIntensity(t,2)=concentration_2(centroid_diff_2(:,1)<= radii_2(:,1),1);

                    YFP_MeanIntensity(t,2)=meanintensity_2(centroid_diff_2(:,1)== min(centroid_diff_2(:,1)),1);
                    YFP_MedianIntensity(t,2)=medianintensity_2(centroid_diff_2(:,1)== min(centroid_diff_2(:,1)),1);
                    YFP_SumIntensity(t,2)=sumintensity_2(centroid_diff_2(:,1)== min(centroid_diff_2(:,1)),1);
                    YFP_ConcIntensity(t,2)=concentration_2(centroid_diff_2(:,1)== min(centroid_diff_2(:,1)),1);


                else
                    YFP_MeanIntensity(t,2)=NaN;
                    YFP_MedianIntensity(t,2)=NaN;
                    YFP_SumIntensity(t,2)=NaN;
                    YFP_ConcIntensity(t,2)=NaN;
                end

                %                 min(centroid_diff_1);

                % % %                                 %% visualization
                % % %                                 centroid_benchmark = figure;
                % % %                                 subplot(1,2,1)
                % % %                                 imshow(CFP_img, []);
                % % %                                 hold on;
                % % %                                 scatter(centroid_tot(r, 1), centroid_tot(r, 2),".");
                % % %                                 hold on;
                % % %                                 scatter(centroid_wc_1(centroid_diff_1==min(centroid_diff_1), 1), centroid_wc_1(centroid_diff_1==min(centroid_diff_1), 2),"red");
                % % %                                 hold off;
                % % %                                 title(sprintf('t = %d CFP Image',t));
                % % %                                 subplot(1,2,2)
                % % %                                 imshow(YFP_img, []);
                % % %                                 hold on;
                % % %                                 scatter(centroid_tot(t, 1), centroid_tot(t, 2),".");
                % % %                                           hold on;
                % % %                                 scatter(centroid_wc_2(:, 1), centroid_wc_2(:, 2),"red");
                % % %                                 hold off;
                % % %                                 title('YFP Image');
                % % %                 % %                 saveas(centroid_benchmark,sprintf('/Users/yijiachen/Documents/B_Cell_project/HVN037/images/benchmark/s1/c%d/centroid_benchmark_c%d_t%d.png',cell_num,cell_num,t))
                % % %                 % %             close all



                % % %                 % Extract centroids from the tracking data and image data
                % % %                 tracked_centroid = centroid_tot(r,:);
                % % %                 CFP_image_centroids = CFP_stats.Centroid;  % From regionprops for CFP
                % % %                 YFP_image_centroids = YFP_stats.Centroid;  % From regionprops for YFP
                % % %
                % % %                 % Plot CFP image with overlaid centroids
                % % %                 figure;
                % % %                 imshow(CFP_img, []);  % Adjust the display range as needed
                % % %                 hold on;
                % % % % % %
                % % % % % % %                 scatter(centroid_tot(:,1), centroid_tot(:,2), 'yellow', 'filled');  % Tracked centroid in red
                % % % % % % %                 hold on
                % % %                 scatter(CFP_image_centroids(:, 1), CFP_image_centroids(:, 2), 'g', 'filled');  % CFP centroids in green
                % % %                 hold on
                % % % % % %                 scatter(centroid_tot([1:10:911],1), centroid_tot([1:10:911],2), 'r', 'x');
                % % %                  scatter(tracked_centroid(1), tracked_centroid(2), 'r', 'x');  % Tracked centroid in red
                % % % % % % %                 hold off;
                % % % % % %                 title(sprintf('CFP Image with Centroids at Frame %d', r));
                % % % % % %                 legend('Tracked Centroid', 'CFP Centroids');
                % % % % % % %                 hold on
                % % % % % %
                % % % % % %
                % % % % % % % % %
                % % % % % % % % % % % %                 % Plot YFP image with overlaid centroids
                % % % % % % % % % % % %                 figure;
                % % % % % % % % % % % %                 imshow(YFP_img, []);  % Adjust the display range as needed
                % % % % % % % % % % % %                 hold on;
                % % % % % % % % % % % %                 scatter(tracked_centroid(1), tracked_centroid(2), 'r', 'filled');  % Tracked centroid in red
                % % % % % % % % % % % %                 scatter(YFP_image_centroids(:, 1), YFP_image_centroids(:, 2), 'g', 'filled');  % YFP centroids in green
                % % % % % % % % % % % %                 hold off;
                % % % % % % % % % % % %                 title(sprintf('YFP Image with Centroids at Frame %d', r));
                % % % % % % % % % % % %                 legend('Tracked Centroid', 'YFP Centroids');


                clear sumintensity_1 sumintensity_2
                clear area_1 area_2
                clear centroid_diff_1 centroid_diff_2
                clear centroid_wc_1 centroid_wc_2
                clear diameters_1 diameters_2
                clear concentration_1 concentration_2
                clear meanintensity_1 meanintensity_2
                clear medianintensity_1 medianintensity_2


            end
        end

        valid_indices_CFP = ~isnan(CFP_MeanIntensity(:, 2)) & CFP_MeanIntensity(:, 2) ~= 0;

        %% mean int
        YFP_CFP_wc_meanint_plot = figure;


        plot(CFP_MeanIntensity(valid_indices_CFP, 1), CFP_MeanIntensity(valid_indices_CFP, 2), 'LineWidth', 2);

        hold on
        valid_indices_YFP = ~isnan(YFP_MeanIntensity(:, 2)) & YFP_MeanIntensity(:, 2) ~= 0;

        plot(YFP_MeanIntensity(valid_indices_YFP, 1), 3.5 * YFP_MeanIntensity(valid_indices_YFP, 2), 'LineWidth', 2);

        ylim([2000 5000])
        xlim([0 960])
        xticks([1 60:60:960])
        % %     %     xticklabels({'0','2','4','6','8','10','12','14','16'})
        xticklabels({'0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16'});

        title('YFP/CFP wholecell mean intensity')
        xlabel('time (h)')
        ylabel('YFP/CFP wholecell mean intensity')
        legend ('CFP','YFP')
        saveas(YFP_CFP_wc_meanint_plot,sprintf('/Users/yijiachen/Documents/B_Cell_project/HVN037/images/fluor_wc_traj/s%d/fluor_plot_meanint_c%d.png',s,cell_num))

        %% median int
        YFP_CFP_wc_medianint_plot = figure;
        valid_indices_CFP = ~isnan(CFP_MedianIntensity(:, 2)) & CFP_MedianIntensity(:, 2) ~= 0;

        plot(CFP_MedianIntensity(valid_indices_CFP, 1), CFP_MedianIntensity(valid_indices_CFP, 2), 'LineWidth', 2);

        hold on
        valid_indices_YFP = ~isnan(YFP_MedianIntensity(:, 2)) & YFP_MedianIntensity(:, 2) ~= 0;

        plot(YFP_MedianIntensity(valid_indices_YFP, 1), 3.5 * YFP_MedianIntensity(valid_indices_YFP, 2), 'LineWidth', 2);

        ylim([2000 5000])
        xlim([0 960])
        xticks([1 60:60:960])
        xticklabels({'0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16'});

        title('YFP/CFP wholecell median intensity')
        xlabel('time (h)')
        ylabel('YFP/CFP wholecell median intensity')
        legend ('CFP','YFP')
        saveas(YFP_CFP_wc_medianint_plot,sprintf('/Users/yijiachen/Documents/B_Cell_project/HVN037/images/fluor_wc_traj/s%d/fluor_plot_medianint_c%d.png',s,cell_num))

        %% sum int
        YFP_CFP_wc_sumint_plot = figure;
        valid_indices_CFP = ~isnan(CFP_SumIntensity(:, 2)) & CFP_SumIntensity(:, 2) ~= 0;

        plot(CFP_SumIntensity(valid_indices_CFP, 1), CFP_SumIntensity(valid_indices_CFP, 2), 'LineWidth', 2);

        hold on
        valid_indices_YFP = ~isnan(YFP_SumIntensity(:, 2)) & YFP_SumIntensity(:, 2) ~= 0;

        plot(YFP_SumIntensity(valid_indices_YFP, 1), 3.5 * YFP_SumIntensity(valid_indices_YFP, 2), 'LineWidth', 2);

        ylim([0.5*10^6 2.5*10^6])
        xlim([0 960])
        xticks([1 60:60:960])
        % %     %     xticklabels({'0','2','4','6','8','10','12','14','16'})
        xticklabels({'0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16'});

        title('YFP/CFP wholecell sum intensity')
        xlabel('time (h)')
        ylabel('YFP/CFP wholecell sum intensity')
        legend ('CFP','YFP')
        saveas(YFP_CFP_wc_sumint_plot,sprintf('/Users/yijiachen/Documents/B_Cell_project/HVN037/images/fluor_wc_traj/s%d/fluor_plot_sumint_c%d.png',s,cell_num))

        %% concentration int
        YFP_CFP_wc_concint_plot = figure;
        valid_indices_CFP = ~isnan(CFP_ConcIntensity(:, 2)) & CFP_ConcIntensity(:, 2) ~= 0;

        plot(CFP_ConcIntensity(valid_indices_CFP, 1), CFP_ConcIntensity(valid_indices_CFP, 2), 'LineWidth', 2);

        hold on
        valid_indices_YFP = ~isnan(YFP_ConcIntensity(:, 2)) & YFP_ConcIntensity(:, 2) ~= 0;

        plot(YFP_ConcIntensity(valid_indices_YFP, 1), 3.5 * YFP_ConcIntensity(valid_indices_YFP, 2), 'LineWidth', 2);

        %                 ylim([0.5*10^6 2.5*10^6])
        xlim([0 960])
        xticks([1 60:60:960])
        % %     %     xticklabels({'0','2','4','6','8','10','12','14','16'})
        xticklabels({'0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16'});

        title('YFP/CFP wholecell concentration intensity')
        xlabel('time (h)')
        ylabel('YFP/CFP wholecell concentration intensity')
        legend ('CFP','YFP')
        saveas(YFP_CFP_wc_concint_plot,sprintf('/Users/yijiachen/Documents/B_Cell_project/HVN037/images/fluor_wc_traj/s%d/fluor_plot_concint_c%d.png',s,cell_num))

        %% save data
        CFP_MedianIntensity_wc_tot{tot_cell_num,1} = CFP_MedianIntensity;
        YFP_MedianIntensity_wc_tot{tot_cell_num,1} = YFP_MedianIntensity;

        CFP_MeanIntensity_wc_tot{tot_cell_num,1} = CFP_MeanIntensity;
        YFP_MeanIntensity_wc_tot{tot_cell_num,1} = YFP_MeanIntensity;

        CFP_SumIntensity_wc_tot{tot_cell_num,1} = CFP_SumIntensity;
        YFP_SumIntensity_wc_tot{tot_cell_num,1} = YFP_SumIntensity;

        CFP_ConcIntensity_wc_tot{tot_cell_num,1} = CFP_ConcIntensity;
        YFP_ConcIntensity_wc_tot{tot_cell_num,1} = YFP_ConcIntensity;

        tot_cell_num = tot_cell_num +1;

        close all
    end

end

save('/Users/yijiachen/Documents/B_Cell_project/HVN037/data_saved/HVN037_CFP_YFP_wc_median_mean_int.mat','CFP_MedianIntensity_wc_tot','YFP_MedianIntensity_wc_tot','CFP_MeanIntensity_wc_tot','YFP_MeanIntensity_wc_tot','CFP_SumIntensity_wc_tot','YFP_SumIntensity_wc_tot','CFP_ConcIntensity_wc_tot','YFP_ConcIntensity_wc_tot')

