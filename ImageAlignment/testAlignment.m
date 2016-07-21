function testAlignment()
%createAlignment creates an alignment mapping out of images to be used with
%worm segmentation and gcamp signal, only uses projective mapping

choice = menu('Image Setup','Split','Multiple');
display('Select registration');
[regFile]=uipickfiles('FilterSpec','Y:\CommunalCode\3dbrain\Registration');
regFile=regFile{1};
registration=load(regFile);


%pick files
display('Select Images');
[fileName]=uipickfiles('FilterSpec','O:\');



if choice==1
    nImage=length(fileName);
else

    nImage=length(fileName)/2;
end

%%
for iImage=1:nImage
    if choice==1
rect1=registration.rect1;
rect2=registration.rect2;
        
    initialIm=double(imread([fileName{iImage}],'tif'));
channelSegment=initialIm((rect1(2)+1):rect1(4),(1+rect1(1)):rect1(3));
channelActivity=initialIm((rect2(2)+1):rect2(4),(1+rect2(1)):rect2(3));
    else
    channelSegment=double(imread([fileName{iImage}],'tif'));
    channelActivity=double(imread([fileName{nImage+iImage}],'tif'));
    end
    
channelSegment=normalizeRange(double(channelSegment(:,:,1)));
channelActivity=normalizeRange(double(channelActivity(:,:,1)));


channelActivity=imwarp(channelActivity,registration.t_concord,'OutputView',...
    registration.Rsegment);

figure
subplot(1,2,1);
imagesc(channelSegment);
subplot(1,2,2);
imagesc(channelActivity);
            
end


