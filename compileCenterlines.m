
function compileCenterlines(dataFolder)

% compileCenterlines works on a dataFolder, find the CLfiles created but
% clusterWormCenterlines, bundles them and calculates the eigenworm
% projections, and centerlines
d= dir([dataFolder filesep 'LowMagBrain*']);
aviFolder=[dataFolder filesep d(1).name];
display(aviFolder)
outputFolder=[aviFolder filesep 'CL_files'];
CLfiles=dir([outputFolder filesep  'CL*']);
CLcell=cell(1,length(CLfiles));
CL_I=cell(1,length(CLfiles));
nCells=16;
for iFile=1:length(CLfiles)
    CLdata=load([outputFolder filesep CLfiles(iFile).name]);
    CLpos=str2double(CLfiles(iFile).name(4:5));
    CLcell{CLpos}=CLdata.CLall;
    CL_I{CLpos}=CLdata.IsAll;
end



%%
CL_Lenght_Fun=@(x)  squeeze(sum(sqrt(sum(diff(x).^2,2))))';

CLcell_2=reshape(CLcell,2,nCells);
CL_I_2=reshape(CL_I,2,nCells);

CLcell_2(2,:)=cellfun(@(x) flipdim(x,3),CLcell_2(2,:),'uniform',0);
CL_I_2(2,:)=cellfun(@(x) flipdim(x,3),CL_I_2(2,:),'uniform',0);

CL_Iall=cell2mat(CL_I_2(:)');

%%
CL_bright_check=mean(CL_Iall>graythresh(CL_Iall));
CL_Iavg=trimmean(CL_Iall(:,CL_bright_check>.4),30,2);

CL_Isum=cellfun(@(x) sum(bsxfun(@minus,x,CL_Iavg).^2), CL_I_2,'uniform',0);
CL_Csum=cellfun(@(x) squeeze(sum(sqrt(sum((diff(x,2)).^2,2)),1))',CLcell_2,'uniform',0);
CL_Isum=cell2mat(CL_Isum);
CL_Csum=cell2mat(CL_Csum);
CL_Isum=CL_Isum+CL_Csum/2;
% CL_Isum(1,:)=smooth(CL_Isum(1,:),300);
% CL_Isum(2,:)=smooth(CL_Isum(2,:),300);
[~,idx]=min(CL_Isum);

CL_length=cellfun(@(x) CL_Lenght_Fun(x),CLcell_2,'uniform',0);
CL_length=cell2mat(CL_length);
% centerline lengths for forward fitting and backward fitting
CL_f=cell2mat(permute(CLcell_2(1,:),[1 3 2]));
CL_b=cell2mat(permute(CLcell_2(2,:),[1 3 2]));
CL_If=cell2mat(CL_I_2(1,:));
CL_Ib=cell2mat(CL_I_2(2,:));
CL_Iout=CL_If;
CL=CL_f;
CL(:,:,idx==2)=CL_b(:,:,idx==2);
CL_Iout(:,idx==2)=CL_Ib(:,idx==2);

NFrames=size(CL_Iout,2);
%% save the data
centerline=CL;
centerline=flip(centerline,2);
save([dataFolder filesep 'centerline_jn3'],'centerline');
centerline=CL_f;
centerline=flip(centerline,2);

save([dataFolder filesep 'centerline_jnf'],'centerline');
centerline=CL_b;
centerline=flip(centerline,2);
save([dataFolder filesep 'centerline_jnb'],'centerline');


%% do eigenprojections and make behavior folder
CLfile=[dataFolder filesep 'centerline_jn3'];
centerlineFile=load(CLfile);
fieldNames=fieldnames(centerlineFile);
clIdx=cellfun(@(x) ~isempty(x), strfind(fieldNames,'line'));

    centerline=centerlineFile.(fieldNames{clIdx(1)});
%create wormcentered coordinate system
wormcentered=FindWormCentered(centerline);
%project onto eigen basis
eigenWormFile='eigenWorms.mat';
load(eigenWormFile);

    
if size(eigbasis,2)~=size(wormcentered,1)
    eigbasis=imresize(eigbasis,[size(eigbasis,1),size(wormcentered,1)]);
    
end
eigenProj=eigbasis*wormcentered;

%save outputs into behaviorAnalysis folder
behaviorFolder=[dataFolder filesep 'BehaviorAnalysis'];
mkdir(behaviorFolder);
offset=0;
save([behaviorFolder filesep 'centerline'] ,'centerline','eigenProj'...
    ,'wormcentered','offset');




