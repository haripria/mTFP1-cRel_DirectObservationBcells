
imageFolder = '/Volumes/GoogleDrive/Shared drives/Signaling Systems Lab/Mark_Xiang/Microscopy/HVN032_BF_ALL_s1_t0001-t3000'; % Replace with the path to your image folder
imageFormat = '*.jpg'; % Replace with the image file format
load('/Users/yijiachen/Documents/B_Cell_project/HVN032/trackmate_data/HVN032_BF_ALL_s1_t0001_t3000_ORG_Concat_trial1_trackdata.mat'); % Load the centroid data

figure;
h_subplot = subplot(1, 1, 1);

frame = 900;
imageFileName = fullfile(imageFolder, sprintf('HVN032_BF_ALL_s1_t0001-t3000_t%04d_ORG.tif', frame));
[currentImage,map] = imread(imageFileName);
imshow(currentImage,map, 'Parent', h_subplot);
hold on;

for c=[1,3:9]
    for sub = 0:50
        frameData(sub+1,:)= alltrackData{c}(alltrackData{c}.Frame == frame-sub, :);
        centroidsX(sub+1,1) = double(frameData.X(sub+1,1)); % Assuming the X-coordinate is in the 3rd column
        centroidsY(sub+1,1) = double(frameData.Y(sub+1,1)); % Assuming the Y-coordinate is in the 4th column

    end

    plot(centroidsX, centroidsY, 'LineWidth', 2); % Overlay a red line connecting centroids
    hold on

end

saveas(h_subplot,sprintf('/Users/yijiachen/Documents/B_Cell_project/HVN032/images/s1/track_snapshot_t%d.png',frame))
