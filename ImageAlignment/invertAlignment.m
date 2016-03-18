function invertAlignment(dataFile,outputName);
% invertAlignment takes a dataFile saved from the createAlignment program
% used on hiRes split data, and inverts the two transformations, switching
% rect1 with rect2, swapping the points, and recreating the transformation.
% 

load(dataFile)
%%
activityPts=Sall;
segmentPts=Aall;

Aall=activityPts;
Sall=segmentPts;
size1=[rect2(3)-rect2(1),rect2(4)-rect2(2)];
rect2_temp=rect2;
rect1_temp=rect1;
rect2=rect1_temp;
rect1=rect2_temp;
t_concord = fitgeotrans(Aall,Sall,'projective');
Rsegment = imref2d(size1);
%% save 
save(['Y:\CommunalCode\3dbrain\registration' outputName],'rect1','rect2','t_concord'...
    ,'Rsegment','padRegion','initialIm','Sall','Aall','fileName')
    

