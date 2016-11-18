#!/bin/sh
#
# File:   runWormCellTracking.sh
# Author: benbratton
#
# Created on May 26, 2015, 11:56:37 AM
# tell the scheduler which version of matlab to use, put it before $PATH so that it takes precedence
export PATH=/usr/local/matlab-R2013a/bin/:$PATH
# This path is for Della
export PATH=/tigress/LICENSED/matlab-R2014b/bin/:$PATH

# parse matlab paths
FILES=$HOME/scripts/shae-pythonSubmissionScripts/*.path
#echo $FILES
for input in $FILES
do
  #echo "Processing $input file..."  
  while read file
  do
    export MATLABPATH="$MATLABPATH:$HOME/scripts/$file"
  done  < "$input"
done

export MATLABPATH="$MATLABPATH;"

# make sure that the matlab path is full/correct
echo $MATLABPATH
echo $1
echo $SGE_TASK_ID
echo $2

# run the job
# if 3rd input is 1, run on with slurm task array ID
   matlab -nosplash -nodesktop -nodisplay -singleCompThread -r "clusterBotCheckCompiler('$1');exit;"
