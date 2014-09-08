folderList=uipickfiles;
tail='conv3d3.mat';
%%
<<<<<<< HEAD
progressbar;
for iFolder=2:length(folderList);
=======
%progressbar;
parpool('local',8,'IdleTimeout', 600);

%%
for iFolder=1:length(folderList);
>>>>>>> origin
    folderName=folderList{iFolder};
    display(['Starting folder: ' folderName]);
    if exist([folderName filesep 'cell_shape_settings_tri.txt'],'file')
        [image_para,flag,cline_para] = paramsText(fullfile(folderName,'cell_shape_settings_tri.txt'));
        
    else
        
        
        [image_para,flag,cline_para] = paramsText(fullfile(folderName,'cell_shape_settings.txt'));
    end
    fileList=dir([folderName filesep '*' tail]);
    fileList={fileList.name};
    parfor iFile=1:length(fileList);
        try
     %   progressbar(iFolder/length(folderList),iFile/length(fileList))
        fileName=[folderName filesep fileList{iFile}];
        data=load(fileName,'cellF','cellcurve','cellF','polyFit');
        cellF=data.cellF;
        cellcurve=data.cellcurve;
        polyFit=data.polyFit;
        mOld=cellF.mFraction;
         mFraction=membraneFractionFile(fileName,image_para,flag);
         cellF.mFraction=mFraction;
        
        
        for iPoly=1:length(polyFit);
             polyFit(iPoly).polyCurve.kg=polymerCurvatureLocalization( cellcurve.kg,polyFit(iPoly).rc );
             polyFit(iPoly).polyCurve.km=polymerCurvatureLocalization( cellcurve.km,polyFit(iPoly).rc );
             polyFit(iPoly).polyCurve.k1=polymerCurvatureLocalization( cellcurve.k1,polyFit(iPoly).rc );
             polyFit(iPoly).polyCurve.k2=polymerCurvatureLocalization( cellcurve.k2,polyFit(iPoly).rc );
            polyFit(iPoly).Fluor=polymerCurvatureLocalization( cellF.Ff,polyFit(iPoly).rc );
        end
                
        keyboard
        
        parsave(fileName,polyFit,'polyFit')
        parsave(fileName,cellF,'cellF');
        catch
        end
    end
end
