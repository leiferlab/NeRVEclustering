function CLfit=eigenWormFitEnergyEigenCLFit(inputImage,CLold,refIdx1,cm,eigbasis)
%function fits an image with a centerline using eigenworms. 
wc=FindWormCentered(CLold);
refIdx2=[refIdx1 99];
ref_point1=cm;
ref_point2=CLold(end,:);

eigVals=eigbasis(:,2:end-1)*wc;
wc_eig=eigVals'*eigbasis(:,2:end-1);
%eigVals=eigVals(1:4);
eigOut=fminsearch(@(x) eigenWormFitEnergy(inputImage,refIdx2,ref_point1,...
    x(end-1:end)',eigbasis,x(1:end-2),1),[eigVals;ref_point2']);

[~,CLfit]=eigenWormFitEnergy(inputImage,refIdx2,ref_point1,...
    eigOut(end-1:end)',eigbasis,eigOut(1:end-2),1);
CLfit=distanceInterp(CLfit,length(CLold));