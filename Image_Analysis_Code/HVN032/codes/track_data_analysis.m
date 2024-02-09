% Specify the file path to your Trackmate data CSV file
addpath '/Users/yijiachen/Documents/B_Cell_project/HVN032/trackmate_data'
filePath_all = '/Users/yijiachen/Documents/B_Cell_project/HVN032/trackmate_data/HVN032_BF_All_s9_t0001-t3000_ORG_concat_spots.csv';  % Replace with the actual file path
filePath_track_result = '/Users/yijiachen/Documents/B_Cell_project/HVN032/trackmate_data/HVN032_BF_All_s9_t0001-t3000_ORG_concat_track_results.xlsx';

% Use readtable to load the data into a table
% Use readtable with 'ReadVariableNames' set to false
trackmateData = readtable(filePath_all);
trackresult= readtable(filePath_track_result);


uniqueCell = unique(trackresult.Cell);
trackData_stitched = cell(length(uniqueCell), 1);

for i = 1:length(uniqueCell)
    cellnum = uniqueCell(i,1);
    rownum = (find(trackresult.Cell == cellnum));
    frame_start = trackresult.frame_start(rownum(1,1));
    frame_end = trackresult.frame_end(rownum(1,1));
    track_ID = trackresult.TrackID(rownum(1,1));

    trackData_stitched = trackmateData(trackmateData.TrackID == track_ID & trackmateData.Frame >= frame_start & trackmateData.Frame <= frame_end,:);
    if size(rownum,1) > 1
        for j = 2:length(rownum)
            frame_start = trackresult.frame_start(rownum(j,1));
            frame_end = trackresult.frame_end(rownum(j,1));
            track_ID = trackresult.TrackID(rownum(j,1));
            trackData_stitched = [trackData_stitched;trackmateData(trackmateData.TrackID == track_ID& trackmateData.Frame >= frame_start & trackmateData.Frame <= frame_end,:)];
        end
    else
    end
    alltrackData{i} = trackData_stitched;
end

save('/Users/yijiachen/Documents/B_Cell_project/HVN032/trackmate_data/HVN032_BF_All_s9_t0001_t3000_ORG_concat_trackdata.mat', 'alltrackData');
