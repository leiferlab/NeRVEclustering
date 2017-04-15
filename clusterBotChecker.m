function clusterBotChecker(filePath,startIdx,nSubSample)
%cluster bot checker uses a path to pointstats and a track index and stepsize as
%inputs. The code goes through all frames and compares them to nSubSample
%number of other randomly selected frames. In each frame, it asks where
%the otherframes would guess where a certain neuron, given by the startIdx,
%should be using TPS interp.

%Inputs:
% filepath : path to pointstats2 file made after clusterTrackCompiler
% startIdx : the number of the run. The program is run for each neuron
% cluster detected in the worm. with <groupSize> number of startIdx devoted
% for each neuron
% groupSize: number of groups to run for each neuron. For example, if a
% video is 2000 volumes long, each index will run 2000/groupSize
% comparison. For 150 neurons, the program will need to be run with
% startIdx up to 300 to do every comparison for every neuron. 

% changed to run under the 1 hr, but no longer needed with della

warning('off');
if nargin==0
    filePath=uipickfiles;
    startIdx=1;
    filePath=filePath{1};
end

if nargin<3
    nSubSample=100;
end
%% setup output paths
outputName=fileparts(filePath);
if isempty(outputName)
    outputName=pwd;
end


outputName=[outputName filesep 'botCheckFolder'];
mkdir(outputName);
%%
%load pointStats file
fileData=load(filePath);
fileField=fieldnames(fileData);
pointStats2=fileData.(fileField{1});

%%
%select number of points,
nTime=length(pointStats2);
%potentially only take subset of times for TPS guessing, for now take all times
subSample=round(1:nTime/nSubSample:nTime);

%which indexes of pointstats to run
startIdxReal=ceil(startIdx);
%which iteration of that pointstats, number between 0 and groupSize-1
runIdx=0;

%list of times
timeVector=1:nTime;


%% loop through selected neurons in list
for iPointIdx=startIdxReal;
    display([' Starting ' num2str(iPointIdx), 'Run' num2str(runIdx)]);
    %initialize comparison matrices. Each of these contains the "guessed"
    %x, y, and z position of each neuron in each frame. 
    comparePointEstimate_x=nan(nSubSample,nTime);
    comparePointEstimate_y=nan(nSubSample,nTime);
    comparePointEstimate_z=nan(nSubSample,nTime);
    comparePointConf=nan(nSubSample,nTime);
    %initialize list of all points
    xyzRefAll=nan(nTime,3);
    
    %loop through all time points
    for iFrame=timeVector
        %%

        %get pointStats being analyzed
        pointsRef=pointStats2(iFrame);
        %get track IDs for each neuron
        refTrackIdx=pointsRef.trackIdx;
        ref_points=pointsRef.straightPoints;
        %populate list of all initially found points
        if any(refTrackIdx==iPointIdx)
            xyzRef=pointsRef.straightPoints(refTrackIdx==iPointIdx,:);
            xyzRefAll(iFrame,:)=xyzRef;
        end
        
        %loop through comparison time points
        for i_samp=1:nSubSample
            %loop through other points and guess where the position of a
            %certain neuron should be
            try
                %get other points and other point track IDs
                pointsI=pointStats2(subSample(i_samp));
                iTrackIdx=pointsI.trackIdx;
                current_points=pointsI.straightPoints;
                %find which points overlap between the current track and
                %the reference track index
                [~,overI,overRef]=intersect(iTrackIdx,refTrackIdx);
                
                %exclude the current point being tested to find the mapping
                %in both current track and ref track
                overI_excludeI=overI(iTrackIdx(overI)~=iPointIdx);
                overRef_excludeI=overRef(refTrackIdx(overRef)~=iPointIdx);
                
                if any(overI) && any(iTrackIdx==iPointIdx)
                    %get points excluding the point in question
                    movingPoints_excludeI=current_points(overI_excludeI,:);
                    controlPoint_excludeI=ref_points(overRef_excludeI,:);
                    xyzI=current_points(iTrackIdx==iPointIdx,:);
                    %also take the confidence of that point
                    xyzI_conf=pointsI.trackWeights(iTrackIdx==iPointIdx);
                    
                    %estimate where the moving frame thinks the reference 
                    %frame's point should be
                    point_estimate=tpswarp3points(...
                        movingPoints_excludeI,controlPoint_excludeI,xyzI);
                    
                    %fill in that point in the matrix of guessed points
                    comparePointEstimate_x(i_samp,iFrame)=point_estimate(1);
                    comparePointEstimate_y(i_samp,iFrame)=point_estimate(2);
                    comparePointEstimate_z(i_samp,iFrame)=point_estimate(3);
                    %also fill in confidence of that point
                    comparePointConf(i_samp,iFrame)=xyzI_conf;
                end
                %display sstatus
                if ~mod(i_samp, 100)
                    display(['Done ' num2str(iFrame), ' and ' num2str(i_samp)])
                end
            catch me
                display(['Error at ' num2str(iFrame), ' and ' num2str(i_samp)])
            end
        end
        
    end
    %save results
    outputNameFile=[outputName filesep 'botChecker' num2str(iPointIdx,'%3.5d')...
        'Run' num2str(runIdx,'%3.2d')];
    display(outputNameFile)
    save(outputNameFile,'comparePointEstimate_x','comparePointEstimate_y', ...
        'comparePointEstimate_z','comparePointConf','xyzRefAll');
    
    
end