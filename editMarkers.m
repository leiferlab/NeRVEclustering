function markerTable = editMarkers(markerFile, centroids)
%EDITMARKERS takes a marker file and creates a new marker file with the 
%(x,y,z) coordinates of each marker changed to those of the center of the
%corresponding nucleus.


%Load the marker file as a table ("markerTable").
markerTable = readtable(markerFile,'FileType','text','TreatAsEmpty',' ');
markerTable.Properties.VariableNames{1} = 'x';


%Extract (x,y,z) coordinates of the markers.
coords = zeros(markerTable.height,3);
coords(:,1) = markerTable.x;
coords(:,2) = markerTable.y;
coords(:,3) = markerTable.z;


%Find the coordinates of the closest nucleus to each marker position.
for i = 1:markerTable.height      %for each marker
    
    %initialize references
    dmin = Inf;
    entry = 0;
    
    for j = 1:centroids.length    %for each centroid
        
        %calculate distance
        x = (coords(i,1) - centroids(j,1))^2;
        y = (coords(i,2) - centroids(j,2))^2;
        z = (coords(i,3) - centroids(j,3))^2;
        d = sqrt(x + y + z);
        
        %update references for likely candidate
        if dmin > d               
            dmin = d;
            entry = j;
        end
    end
    coords(i,:) = centroids(entry,:);
end


%Replace the coordinates of the marker with those of the centrer of the
%closest nucleus in the table.
markerTable.x = coords(:,1);
markerTable.y = coords(:,2);
markerTable.z = coords(:,3);


%Create a new marker file with the corrected coordinates
oldFile = markerFile(1:(length(markerFile)-7));
newFile = strcat(oldFile,'_centered.csv');

writetable(markerTable, newFile);


end