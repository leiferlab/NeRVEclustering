#!/bin/sh
#
# File:   runMatlabInput.sh
# Author: jnguyen

# runs any input string in matlab, designed for use with cluster
#
# Created on May 26, 2015, 11:56:37 AM
# tell the scheduler which version of matlab to use, put it before $PATH so that it takes precedence
# this path for genomics
export PATH=/usr/local/MATLAB/R2015a/bin/:$PATH
# This path is for Della
export PATH=/tigress/LICENSED/matlab-R2017a/bin/:$PATH

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

# export path
export MATLABPATH="$MATLABPATH;"

# make sure that the matlab path, input are full/correct
echo $MATLABPATH
echo $1
# run the job, supressing figures
matlab -nosplash -nodesktop -nodisplay -singleCompThread -r "$1;exit;"
