#!/bin/sh
#
# File:   runWormCenterlineFitting.sh
# Author: Jeffrey Nguyen
#
# Created on May 26, 2015, 11:56:37 AM
# tell the scheduler which version of matlab to use, put it before $PATH so that it takes precedence
export PATH=/usr/local/matlab-R2013a/bin/:$PATH
# This path is for Della
export PATH=/tigress/LICENSED/matlab-R2014b/bin/:$PATH

# parse matlab paths
FILES=$CODE_HOME/3dbrain/PythonSubmissionScripts/*.path
#echo $FILES
for input in $FILES
do
  #echo "Processing $input file..."  
  while read file
  do
    export MATLABPATH="$MATLABPATH:$CODE_HOME/$file"
  done  < "$input"
done

export MATLABPATH="$MATLABPATH;"

# make sure that the matlab path is full/correct
echo $MATLABPATH
# run the job
   matlab -nosplash -nodesktop -nodisplay -singleCompThread -r "clusterWormCenterline('$1',$SLURM_ARRAY_TASK_ID);exit;"
