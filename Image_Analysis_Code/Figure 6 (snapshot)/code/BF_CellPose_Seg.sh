


#!/bin/bash
# for the following CellPose Pipeline, please run in the virtual environment.
#conda activate cellpose

# Set the name of the log file and the path where it should be saved
LOG_FILE=/Users/yijiachen/Documents/B_Cell_project/panel_cd/data_saved/logfile.log

# Temporary directory to store CellPose output
TEMP_DIR=/Users/yijiachen/Documents/B_Cell_project/panel_cd/temp
IMG_DIR=/Users/yijiachen/Documents/B_Cell_project/panel_cd/HVN032_Fluo_9_scenes/HVN032_Fluo_ALL_s$(printf "%01d" $s)
#set up the scene number s
for s in $(seq 1 1 3)
do
for i in $(seq 1 2 150)
#1 2 3 4 7 13 19 25 49 97 145
#$(seq 154 1 3000)
do
    in_file=/Users/yijiachen/Documents/B_Cell_project/panel_cd/HVN032_Fluo_9_scenes/HVN032_Fluo_ALL_s$(printf "%01d" $s)/HVN032_Fluo_ALL_t$(printf "%03d" $i)c1_ORG.tif
    echo $in_file
    diameter=$(($(($i/500))+18))
    echo $diameter
    # Run the Python script and store output in the temporary directory
    /Users/yijiachen/anaconda3/bin/python -m cellpose --verbose --use_gpu --image_path $in_file  --pretrained_model cyto2 --net_avg --diameter $diameter --min_size 15 --save_png --savedir $TEMP_DIR >> $LOG_FILE 2>&1
    
    # Move only the mask image to the final destination and rename it
    mv $TEMP_DIR/HVN032_Fluo_ALL_t$(printf "%03d" $i)c1_ORG_cp_masks.png /Users/yijiachen/Documents/B_Cell_project/panel_cd/HVN032_BF_mask/s$(printf "%01d" $s)/HVN032_Fluo_ALL_t$(printf "%03d" $i)_cp_masks.png
   
    # Clean up the temporary directory
    rm -f $TEMP_DIR/HVN032_Fluo_ALL_t$(printf "%03d" $i)_cp_output.png
    rm -f $IMG_DIR/HVN032_Fluo_ALL_t$(printf "%03d" $i)c1_ORG_seg.npy
done
done
# Remove the temporary directory when done
rmdir $TEMP_DIR


