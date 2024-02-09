#!/bin/bash

for scene in 2 3 4 5 6 7 8
do
    maskdir=/home/mark/cellpose/HVN037a_BF_ALL_cellpose_masks
    mkdir -p $maskdir/s$scene\_merged

    # Define the path to the ImageJ executable
    imagej="./ImageJ-linux64"

    # Define the path to the macro script
    macrodir=/home/mark/cellpose

    LOG_FILE_Merge=$maskdir/logfile.log

    # # Run ImageJ in headless mode with the macro script for each frame
    # for i in $(seq 1 1 2)
    # do
    #     macro_arg=$scene+$i
    #     mask_file=$maskdir/s$scene\_merged/HVN032_BF_ALL_s$scene\_t0001-t3000_t$(printf "%04d" $i)_ORG_Merged.tif
    #     if [ ! -f $mask_file ]; then
    #         #run macro with scene and frame input
    #         $imagej --headless --console -macro $macrodir/Merge_sb.ijm $macro_arg
    #     fi

    #     # if [ ! -f $mask_file ]; then
    #     #     echo "Failed to add: $mask_file" >> $LOG_FILE_Merge 2>&1
    #     #     #run macro with scene and frame input
    #     #     $imagej --headless --console -macro $macrodir/Merge_sb.ijm $macro_arg
    #     # fi
    # done

    $imagej --headless --console -macro $macrodir/Concat_sb.ijm $scene
done

