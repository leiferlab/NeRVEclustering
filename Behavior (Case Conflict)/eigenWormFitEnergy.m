

function [E,centerline_eig]=eigenWormFitEnergy(inputImage,refIdx2,ref_point1,ref_point2,eigbasis,eigVals,show)
if nargin==5
    show=0;
end


wc_eig=eigVals'*eigbasis(1:length(eigVals),2:end-1);
centerline_eig=recreateWorm(wc_eig,refIdx2,[ref_point1;ref_point2],500);
CLX=centerline_eig(:,1);
CLY=centerline_eig(:,2);
E=CLsearch(inputImage,CLX,CLY,show);


