function angle=GetRotateMask(im,FourierMask)
% calculate the angle to rotate the fourier mask.

%% rotation
    %figure,imagesc(im);
    fftimg=fftshift(fft2(im));
    original=abs(fftimg);
    %figure,imagesc(original)
    stdfilimg=stdfilt(original);
    %figure,imagesc(stdfilimg);
    
    angstep=40;
    theta=(0:angstep)/angstep*180;
    correlation=zeros(angstep+1,1);
    for i=1:(angstep+1)
        rotateimg=imrotate(stdfilimg,theta(i),'crop');
        %figure,imagesc(rotateimg);
        correlation(i)=corr2(rotateimg,FourierMask)+1.147e-4*(i-1); %correction for image rotate.
    end
    
    [~,maxcoridx]=max(correlation);
    angle=-theta(maxcoridx);