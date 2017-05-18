function Demo_pathsetup()

%simple code made for adding scripts to the path. This is necessary for
%running any of the democode on the sample datasets released with Nguyen et
%al 2017. 


file_name=mfilename('fullpath');
p=fileparts(fileparts(file_name));
addpath(genpath(p));
