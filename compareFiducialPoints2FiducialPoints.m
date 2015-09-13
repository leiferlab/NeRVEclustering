%%%% compare two sets of fiducial point cells by finding correspondences in
%%%% all frames
display('Select both fiducial point files')
fFiles=uipickfiles();
%%
fiducials1=load(fFiles{1});
fiducials2=load(fFiles{2});
fiducials1=fiducials1.fiducialPoints;
fiducials2=fiducials2.fiducialPoints;


%% tracking params
param.dim=3;
param.excessive=4;
 param.quiet=1;
 param.difficult=5.e3;
 param.good=2;
 
 minDist=10;

%%
searchList=find(cellfun(@(x) ~isempty(x),fiducials1));

matchMatrix=nan(150,max(searchList));

for i=searchList
    display(['starting round' num2str(i)])
    points1=fiducials1{i};
    points2=fiducials2{i};
    if ~isempty(points1)
        empty1=cellfun(@(x) isempty(x),points1);
        empty1list=find(all(~empty1,2));
        points1(empty1)={nan};
        empty2=cellfun(@(x) isempty(x),points2);
        empty2list=find(all(~empty2,2));
        points2(empty2)={nan};  
        points1=cell2mat(points1);
        points2=cell2mat(points2);
        trackIn1=[points1(empty1list,:) empty1list ones(size(empty1list))];
        trackIn2=[points2(empty2list,1:end-1) empty2list 2*ones(size(empty2list))];
        
        trackInput=[trackIn1;trackIn2];
        trackInput=trackInput(:,[1 2 4 5 end]);
        trackInput(:,3)=trackInput(:,3)-min(trackInput(:,3))+1;
        trackInput(:,3)=trackInput(:,3)*10;
        trackOutput=nan;
        trackCounter=0;
        while isnan(trackOutput)
        trackOutput=trackJN(trackInput,20-trackCounter,param);
        trackCounter=trackCounter+1
        
        end
        paired1=trackOutput(1:2:end,4);
        paired2=trackOutput(2:2:end,4);
        matchMatrix(paired1,i)=paired2;
    end
    
    
end

%%
matchMatrixRev=nan(size(matchMatrix));

for i=1:size(matchMatrix,2);
    paired1=find(~isnan(matchMatrix(:,i)));
        if ~isempty(paired1)

    paired2=matchMatrix(paired1,i);
   matchMatrixRev(paired2,i)=paired1;
    end
    
end

%%

currentMatrix=matchMatrix;
corrMatch=bsxfun(@eq,currentMatrix,mode(currentMatrix,2));
wrongMatch=~corrMatch & ~isnan(currentMatrix);
viewMatch=corrMatch*2+wrongMatch;
viewMatch=sort(viewMatch,2);
