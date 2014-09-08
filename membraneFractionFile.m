function mFraction=membraneFractionFile(fileName,image_para,flag)
load(fileName,'-regexp','^(?!eg$)...');

if strfind(fileName, 'TRI')
    triflag=1;
else
    triflag=0;
end

if triflag
[V1,V2]=stackLoad([fileName(1:strfind(fileName,'TRI')-1) '.tif'],[],image_para.ZX_factor);
else
    [V1,V2]=stackLoad([fileName(1:strfind(fileName,'conv3d3')-1) '.tif'],[],image_para.ZX_factor);

end

if flag.reverse
    hyper_stack=V2;
    F1=V1;
else
    hyper_stack=V1;
    F1=V2;
end
        image_para.imsize=size(F1);




F1=normalizeRange(pedistalSubtract(F1));
hyper_stack=normalizeRange(pedistalSubtract(hyper_stack));

if isfield(flag,'F1fix')
    if flag.F1fix
        AAA=convnfft(hyper_stack,...
            F1(end:-1:1,end:-1:1,end:-1:1));
        [xfix,yfix,zfix]=ind2sub(size(AAA),find(AAA==max(AAA(:))));
        xfix=xfix-ceil(size(AAA,1)/2);
        yfix=yfix-ceil(size(AAA,2)/2);
        zfix=zfix-ceil(size(AAA,3)/2);
        
        
        [Xmap,Ymap,Zmap]=ndgrid(1:size(F1,1),1:size(F1,2),1:size(F1,3));
        Xmap=Xmap-xfix;
        Ymap=Ymap-yfix;
        Zmap=Zmap-zfix;
        F1=interp3(F1,Ymap,Xmap,Zmap,'nearest',0);
    end
end




if triflag
    tri=F;
    
    
    averageSurfPSF = getappdata(0,'averageSurfPSF');
    aPSF=normalizeRange(averageSurfPSF(:,:,:,1));
    aPSF=image_resize(aPSF,size(aPSF,1),...
        size(aPSF,2),round(image_para.ZX_factor*size(aPSF,3)));
    
    tri=triInterp2(tri,[],2);
    membrane2=coord2image3d(tri,image_para.imsize,1,image_para.nm_per_pixel);
    mfill=imfill(double(membrane2~=0),'holes');
    convMem=hyper_stack;%convnfft(membrane2,aPSF,'same');
    convMem=normalizeRange(convMem);
    convFill=convnfft(mfill,aPSF,'same');
    convFill=normalizeRange(convFill);
    
    x=fminsearch(@(x) sum(sum(sum((F1-(x(1)*convFill+(x(2))*convMem)).^2))),[.5,.5]);
    convFill=convFill.*x(1);
    Fmembrane=F1-convFill;
    Fmembrane(Fmembrane<0)=0;
    Ffill=max(sum(convFill(:)),0);
    Fmem=max(sum(Fmembrane(:)),0);
    mFraction=Fmem/(Fmem+Ffill);
    
else
    

    
    
    averageSurfPSF = getappdata(0,'averageSurfPSF');
    
    membrane2=coord2image3d(cellcoord.x,cellcoord.y,cellcoord.z,image_para.imsize,1,image_para.nm_per_pixel);
    mfill=imfill(double(membrane2~=0),'holes');
    close all
    convMem=hyper_stack;%convnfft(membrane2,aPSF,'same');
    convMem=normalizeRange(convMem);
    convFill=convnfft(mfill,averageSurfPSF,'same');
    convFill=normalizeRange(convFill);
    
    [x,res]=fminsearch(@(x) sum(sum(sum((F1-(x(1)*convFill+(x(2))*convMem)).^2))),[.5,.5]);
    convFill=convFill.*x(1);
    %Fmembrane=F1-convFill;
    %Fmembrane(Fmembrane<0)=0;
    convFill=F1-Fmembrane;
    convFill(convFill<0)=0;
    Ffill=max(sum(convFill(:)),0);
    Fmem=max(sum(Fmembrane(:)),0);
    mFraction=Fmem/(Fmem+Ffill);
end









%%

