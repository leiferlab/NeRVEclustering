function FourierMask=GenerateFourierMask(behavior_vidobj)
% use the low mag image  and generate the fourier mask.

nframes=round(behavior_vidobj.Duration*behavior_vidobj.FrameRate);
skip=max(20,round(nframes/100));

rFrames=1:skip:nframes;
theta=zeros(length(rFrames),1);
count=1;
for itime=rFrames
    %progressbar(itime/nframes);
    behavior_frame = read(behavior_vidobj,itime);
    %Default fourier mask without rotation
    if  ~exist('FourierMask','var')
        hexlength=97;
        FourierMask=FourierMaskHexagon(hexlength,behavior_frame,4);
    end
    theta(count)=GetRotateMask(behavior_frame,FourierMask);
    count=count+1;
        
end

angletable=tabulate(theta);
[percent,angle_idx]=max(angletable(:,3));
if percent>60
    FourierMask=imrotate(FourierMask,angletable(angle_idx,1),'crop');
else 
    FourierMask=[];
end


