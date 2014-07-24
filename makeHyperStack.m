function makeHyperStack
%% get a list of folders
folderList = uipickfiles();
%% make the assumption that everything can be read in to memory
% delete('temp2.tif');

for iFolder = 1:size(folderList,2)
    fileList = dir([folderList{iFolder},filesep,'*.tif']);
    for jFile = 1:numel(fileList)
        imageStack = imread(fullfile(folderList{iFolder},fileList(jFile).name));
        if and(jFile==1,iFolder==1);
            imwrite(imageStack,'temp.tif','tif','Compression','none','writeMode','overwrite');
            
        else
            
            imwrite(imageStack,'temp.tif','tif','Compression','none','writeMode','append');
        end
    end
end

disp(['done ',datestr(now)]);
%%
imwrite(0,'temp.tif','tif','Compression','none','writeMode','overwrite');