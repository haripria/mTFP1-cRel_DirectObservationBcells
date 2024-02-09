#!/bin/bash
# for the following CellPose Pipeline, please run in the virtual environment.
#conda activate cellpose

# Set the name of the log file and the path where it should be saved
scene=2
filedir=/home/mark/cellpose/HVN032_BF_ALL_cellpose_masks
outdir=$filedir/s$scene
mkdir -p $outdir

indir=/home/mark/Signaling_Systems_Lab_Google_Drive/Mark_Xiang/Microscopy/HVN032_BF_ALL_s$scene\_t0001-t3000
LOG_FILE=$outdir/logfile.log

for i in $(seq 2296 1 3000)
#for i in 2295
do
    in_file=$indir/HVN032_BF_ALL_s$scene\_t0001-t3000_t$(printf "%04d" $i)_ORG.tif
    #echo $in_file
    # Run the Python script and append the output to the log file
    diameter=$(($(($i/500))+18))
    #echo $diameter
    python -m cellpose --verbose --use_gpu --image_path $in_file  --pretrained_model cyto2 --net_avg --diameter $diameter \
           --min_size 15 --save_png --savedir $outdir >> $LOG_FILE 2>&1
done

