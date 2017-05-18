#!/bin/sh
#
# File:   runWormBotChecker.sh
# Author: jnguyen
#
# Runs the error correction step fro the NeRVE clustering algorithm for tracking neruons.
# Created on May 26, 2015, 11:56:37 AM
# tell the scheduler which version of matlab to use, put it before $PATH so that it takes precedence
export PATH=/usr/local/matlab-R2013a/bin/:$PATH
# This path is for Della
export PATH=/tigress/LICENSED/matlab-R2017a/bin/:$PATH

# parse matlab path files
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
echo $1
echo $2

# run the job
   matlab -nosplash -nodesktop -nodisplay -singleCompThread -r "clusterBotChecker('$1',$SLURM_ARRAY_TASK_ID, $2);exit;"
