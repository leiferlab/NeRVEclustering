function [rows,cols]=getdatdimensions(string)

%getdatdimension takes a string produced by the whole brain imaging system,
%which has the image dimensions in the text, and parses the text to extract
%the rows and columns. string is in the format sCMOS_Frames_U16_1024x512.dat
%  [rows,cols]=getdatdimensions(string)

xloc= find(string=='x',1,'last');
dotloc=find(string=='.',1,'last');
numstart= find(string=='_',1,'last');
rows=str2double(string(numstart+1:xloc-1));
cols=str2double(string(xloc+1:dotloc-1));

