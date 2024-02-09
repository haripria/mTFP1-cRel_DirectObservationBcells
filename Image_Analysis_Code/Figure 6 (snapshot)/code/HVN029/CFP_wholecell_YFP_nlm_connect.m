% Assuming you have loaded both CFP_wholecell_all_t and YFP_nucleus_live_mask_all_t structures
addpath '/Users/yijiachen/Documents/B_Cell_project/panel_cd/data_saved'
load('HVN029_YFP_nucleus_live_mask_value.mat')
load('HVN029_CFP_wholecell_value.mat')


% Loop through the rows of CFP_wholecell_all_t.sX.wc
for s = 1:8
    t_tot = [1 2 3 4 7 13 19 25 49 97 145];
    fn=['s', num2str(s)];
    for time_wc_nlm = 1:size(CFP_wholecell_all_t.(fn).('wc_sumint'), 2)
        row_connected = 1;
       for cell_wc = 1:size(CFP_wholecell_all_t.(fn).('wc_sumint'), 1)
            col_time_point_1 = (2*(time_wc_nlm-1)+1);
            col_time_point_2 = (2*(time_wc_nlm-1)+1)+1;
            % Get the centroid of the current whole cell in CFP_wholecell_all_t
            centroid_wc = CFP_wholecell_all_t.(fn).('wc_centroid')(cell_wc, col_time_point_1:col_time_point_2);

            % Loop through the rows of YFP_nucleus_live_mask_all_t.sX.nucleus (assuming a similar structure exists)

            for cell_nucleus = 1:size(YFP_nucleus_live_mask_all_t.(fn).('nlm_sumint'), 1)
                % Get the centroid of the current nucleus in YFP_nucleus_live_mask_all_t
                centroid_nucleus = YFP_nucleus_live_mask_all_t.(fn).('nlm_centroid')(cell_nucleus, col_time_point_1:col_time_point_2);
                centroidThreshold = YFP_nucleus_live_mask_all_t.(fn).('radius')(cell_nucleus, time_wc_nlm);
                % Calculate the distance between centroids
                centroidDistance = ((centroid_wc(1,1)-centroid_nucleus(1,1))^2 + (centroid_wc(1,2)-centroid_nucleus(1,2))^2)^(1/2);


                % Check if the centroids are within the threshold distance
                if centroidDistance <= centroidThreshold
                    % Pair the intensity values (you can modify this based on your data structure)
                    %% sum intensity
                    wc_sumint = CFP_wholecell_all_t.(fn).('wc_sumint')(cell_wc, time_wc_nlm);
                    nucleus_sumint = YFP_nucleus_live_mask_all_t.(fn).('nlm_sumint')(cell_nucleus, time_wc_nlm);
                    int_connected_CFP_YFP.(fn).('connectedCFP_YFP_sumint_wc')(row_connected,time_wc_nlm)=wc_sumint;
                    int_connected_CFP_YFP.(fn).('connectedCFP_YFP_sumint_nucleus')(row_connected,time_wc_nlm)=nucleus_sumint;
                    %% median intensity
                    wc_medianint = CFP_wholecell_all_t.(fn).('wc_medianint')(cell_wc, time_wc_nlm);
                    nucleus_medianint = YFP_nucleus_live_mask_all_t.(fn).('nlm_medianint')(cell_nucleus, time_wc_nlm);
                    int_connected_CFP_YFP.(fn).('connectedCFP_YFP_medianint_wc')(row_connected,time_wc_nlm)=wc_medianint;
                    int_connected_CFP_YFP.(fn).('connectedCFP_YFP_medianint_nucleus')(row_connected,time_wc_nlm)=nucleus_medianint;
                    %% centroid
                    int_connected_CFP_YFP.(fn).('connectedCFP_YFP_centroid_nucleus')(row_connected,2*(time_wc_nlm-1)+1:2*(time_wc_nlm-1)+2)=centroid_nucleus;
                    int_connected_CFP_YFP.(fn).('connectedCFP_YFP_centroid_wc')(row_connected,2*(time_wc_nlm-1)+1:2*(time_wc_nlm-1)+2)=centroid_wc;

                    row_connected=row_connected+1;
                    % You now have matching intensity values for this pair of cells
                    % You can store or process them as needed
                end
            end
        end
%         %% for visualization
%         t=t_tot(1,time_wc_nlm);
%         path = sprintf('/Users/yijiachen/Documents/B_Cell_project/panel_cd/HVN029_live_nucleus_mask/YFP_masked/s%d', s);
%         addpath(path);
%         filename_YFP_nlm = ['YFP_img_masked_', num2str(t), '.tiff'];
%         info = imfinfo(filename_YFP_nlm);
% 
%         YFP_img_nlm = zeros(info(1).Height,info(1).Width,length(info),'uint16');
%         for read = 1:length(info)
%             YFP_img_nlm(:,:,read) = imread(filename_YFP_nlm,read);
%         end
%         figure
%         imagesc(YFP_img_nlm)
%         truesize([500 500]);
%         colormap gray
%         colorbar
%         title("YFP_img_nlm")
%         hold on
%         scatter(int_connected_CFP_YFP.(fn).('connectedCFP_YFP_centroid_wc')(:,2*(time_wc_nlm-1)+1), int_connected_CFP_YFP.(fn).('connectedCFP_YFP_centroid_wc')(:,2*(time_wc_nlm-1)+2), 'r.', 'SizeData', 100);
%         scatter(int_connected_CFP_YFP.(fn).('connectedCFP_YFP_centroid_nucleus')(:,2*(time_wc_nlm-1)+1), int_connected_CFP_YFP.(fn).('connectedCFP_YFP_centroid_nucleus')(:,2*(time_wc_nlm-1)+2), 'b.', 'SizeData', 50);
%         my_data=int_connected_CFP_YFP.(fn).('connectedCFP_YFP_centroid_nucleus')(:,2*(time_wc_nlm-1)+1:2*(time_wc_nlm-1)+2);
%         % Plot the centroids and row numbers
%         for i = 1:size(my_data,1) % Replace 'your_data' with the appropriate variable containing centroids
%             x = my_data(i, 1); % X-coordinate of the centroid
%             y = my_data(i, 2); % Y-coordinate of the centroid
%             rowNumber = i; % Row number (can be any identifier you want)
% 
%             % Take the integer part of the mean intensity value
%             rowNumber = floor(rowNumber);
% 
%             % Add the integer part of the row number as text at the centroid
%             text(x, y, num2str(rowNumber), 'Color', 'white', 'FontSize', 10, 'FontWeight', 'bold');
%         end
% 
%         %% visualization wholecell CFP
%         path = sprintf('/Users/yijiachen/Documents/B_Cell_project/panel_cd/HVN029_live_nucleus_mask/CFP_wholecell_masked/s%d', s);
%         addpath(path);
%         filename_CFP_wc = ['CFP_wholecell_masked_', num2str(t), '.tiff'];
%         info = imfinfo(filename_CFP_wc);
% 
%         CFP_img_wc = zeros(info(1).Height,info(1).Width,length(info),'uint16');
%         for read = 1:length(info)
%             CFP_img_wc(:,:,read) = imread(filename_CFP_wc,read);
%         end
%         figure
%         imagesc(CFP_img_wc)
%         truesize([500 500]);
%         colormap gray
%         colorbar
%         title("CFP_img_wc")
%         hold on
%         scatter(int_connected_CFP_YFP.(fn).('connectedCFP_YFP_centroid_wc')(:,2*(time_wc_nlm-1)+1), int_connected_CFP_YFP.(fn).('connectedCFP_YFP_centroid_wc')(:,2*(time_wc_nlm-1)+2), 'r.', 'SizeData', 100);
%         scatter(int_connected_CFP_YFP.(fn).('connectedCFP_YFP_centroid_nucleus')(:,2*(time_wc_nlm-1)+1), int_connected_CFP_YFP.(fn).('connectedCFP_YFP_centroid_nucleus')(:,2*(time_wc_nlm-1)+2), 'b.', 'SizeData', 50);
%         my_data=int_connected_CFP_YFP.(fn).('connectedCFP_YFP_centroid_nucleus')(:,2*(time_wc_nlm-1)+1:2*(time_wc_nlm-1)+2);
%         % Plot the centroids and row numbers
%         for i = 1:size(my_data,1) % Replace 'your_data' with the appropriate variable containing centroids
%             x = my_data(i, 1); % X-coordinate of the centroid
%             y = my_data(i, 2); % Y-coordinate of the centroid
%             rowNumber = i; % Row number (can be any identifier you want)
% 
%             % Take the integer part of the mean intensity value
%             rowNumber = floor(rowNumber);
% 
%             % Add the integer part of the row number as text at the centroid
%             text(x, y, num2str(rowNumber), 'Color', 'white', 'FontSize', 10, 'FontWeight', 'bold');
%         end

    end
end

save('/Users/yijiachen/Documents/B_Cell_project/panel_cd/data_saved/HVN029_CFP_wholecell_YFP_nlm_connected.mat','int_connected_CFP_YFP')


