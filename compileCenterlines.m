
function compileCenterlines(dataFolder)

% compileCenterlines works on a dataFolder, find the CLfiles created by
% clusterWormCenterlines, bundles them and calculates the eigenworm
% projections, and centerlines. The centerlines are produced with a forward
% pass and a back pass. The program attempts to pick out the best one. The
% input is the dataFolder being analyzed that contains the dat file. 


%% get folders that have the avi files
%Lowmag folder only exist for new set
d= dir([dataFolder filesep 'LowMagBrain*']);
if isempty(d)
    %for old
    aviFolder=dataFolder;
else
    %for new
    aviFolder=[dataFolder filesep d(1).name];
end

display(aviFolder)
%% load eigenbasis
eigenWormFile='eigenWorms.mat';
load(eigenWormFile);


%% load CL files
outputFolder=[aviFolder filesep 'CL_files'];
CLfiles=dir([outputFolder filesep  'CL*']);
nCells=16;

CLcell=cell(1,2*nCells);
CL_I=cell(1,2*nCells);
for iFile=1:length(CLfiles)
    CLdata=load([outputFolder filesep CLfiles(iFile).name]);
    CLpos=str2double(CLfiles(iFile).name(4:5));
    %load all CL coordinates from each job
    CLcell{CLpos}=CLdata.cl_all;
    %load all CL pixel intensities from each job
    CL_I{CLpos}=CLdata.cl_intensities;
end


%%
%make simple length function 
CL_Lenghth_Fun=@(x)  squeeze(sum(sqrt(sum(diff(x,2).^2,2))))';

%reorganize data to 2x16 cell array, with 1x16 forward, 1x16 backward
CLcell_2=reshape(CLcell,2,nCells);
CL_I_2=reshape(CL_I,2,nCells);
%permute matrix dimensions for ease.
CLcell_2(2,:)=cellfun(@(x) flipdim(x,3),CLcell_2(2,:),'uniform',0);
CL_I_2(2,:)=cellfun(@(x) flipdim(x,3),CL_I_2(2,:),'uniform',0);

%fill in missing jobs with complete job for that cell (replace forward with
%back and vice versa as needed)
isMissing=cellfun(@(x) isempty(x), CL_I_2);
cell_nframes=cellfun(@(x)  size(x,2),CL_I_2);
isMissing= or(isMissing, bsxfun(@ne,cell_nframes,max(cell_nframes,[],2)));
%also replace data that is short for some reason,, will find out exactlyl
%why later but probably means something went wrong


replace=circshift(isMissing,[1,0]);
CLcell_2(isMissing)=CLcell_2(replace);
CL_I_2(isMissing)=CL_I_2(replace);

%convert intensities to matrix.
CL_Iall=cell2mat(CL_I_2(:)');

%% try to determine some measure of how good the CL works
%check intensities
%get mean intensities of pixels above some threshold
CL_bright_check=mean(CL_Iall>graythresh(CL_Iall));
CL_Iavg=trimmean(CL_Iall(:,CL_bright_check>.4),30,2);

%get mean subtracted intensities
CL_Isum=cellfun(@(x) sum(bsxfun(@minus,x,CL_Iavg).^2), CL_I_2,'uniform',0);
%get lengths of each CL
CL_length=cellfun(@(x) CL_Lenghth_Fun(x),CLcell_2,'uniform',0);
%weight for worms that arent too long, and arent too bright (might be
%refelction, subjecto to change)
CL_Isum=cell2mat(CL_Isum);
CL_Csum=cell2mat(CL_length);
CL_Isum=CL_Isum+CL_Csum/100;

% see how curvy the worm is, better worms normally fit the first modes
% beter(regularyly curved)

%reshape eigenbasis for matrix multiplication
eigbasis=imresize(eigbasis,[size(eigbasis,1),size(CLcell_2{1},1)-1]);
wc=cellfun(@(x) sum(eigbasis(1:2,:)*FindWormCentered(x)).^2, CLcell_2,'uniform',0);
wc2=cell2mat(wc);

CL_Isum=CL_Isum-wc2;

%% add something to downweight stuck frames

for iCell=1:size(CLcell_2,2)
CL1=CLcell_2{1,iCell};
CL2=CLcell_2{2,iCell};
I1=(squeeze(mean(sqrt(sum(diff(CL1,[],3).^2,2)))));
I2=(squeeze(mean(sqrt(sum(diff(CL2,[],3).^2,2)))));
I1(I1<0.01)=max(I1);
I2(I2<0.01)=max(I2);
L1=squeeze(sum(sqrt(sum(diff(CL1,[],1).^2,2))));
L2=squeeze(sum(sqrt(sum(diff(CL2,[],1).^2,2))));

I1=[0;cumsum(I1)];
I2=[cumsum(I2,'reverse');0];
I3=squeeze(mean(sqrt(sum((CL1-CL2).^2,2))));

[~,loc]=min(I1+I2+I3);
CL_current=CL1;
CL_current(:,:,loc:end)=CL2(:,:,loc:end);
CL_out{iCell}=CL_current;
end


%% pick out centerlines that look better

%organise matrices
CL_f=cell2mat(permute(CLcell_2(1,:),[1 3 2]));
CL_b=cell2mat(permute(CLcell_2(2,:),[1 3 2]));
CL_If=cell2mat(CL_I_2(1,:));
CL_Ib=cell2mat(CL_I_2(2,:));
CL_Iout=CL_If;

%take only parts of CL's we like
[~,idx]=min(CL_Isum);
CL=cell2mat(permute(CL_out,[1 3 2]));
CL_Iout(:,idx==2)=CL_Ib(:,idx==2);

%% save the data
centerline=CL;
centerline=flip(centerline,2);
save([dataFolder filesep 'centerline_jn3'],'centerline');
%save forward CL
centerline=CL_f;
centerline=flip(centerline,2);
save([dataFolder filesep 'centerline_jnf'],'centerline');
%save reverse CL
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

    

eigenProj=eigbasis*wormcentered;
%save outputs into behaviorAnalysis folder
behaviorFolder=[dataFolder filesep 'BehaviorAnalysis'];
mkdir(behaviorFolder);
save([behaviorFolder filesep 'centerline'] ,'centerline','eigenProj'...
    ,'wormcentered');

%% write to status
hostname = char( getHostName( java.net.InetAddress.getLocalHost ) );
if contains(hostname,'della')
    Fid=fopen([dataFolder filesep 'status.txt'],'a');
    status=[datestr(datetime('now')) ':Finished Centerlines \n'];
    fprintf(Fid,status);
    fclose(Fid);
end

