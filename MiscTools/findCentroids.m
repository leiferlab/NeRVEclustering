function [centroids, numNeurons] = findCentroids(binCoords)
%FINDCENTROIDS takes a 3D matrix representation of a binary image and
%returns an Nx3 matrix of the coordinates of the connected components in
%the image ("centroids").


%Compiles the connected components ("CC") from the binary image
%representation ("binCoords") and finds how many there are ("numNeurons").
CC = bwconncomp(binCoords);
numNeurons = CC.NumObjects;

%Finds the coordinates of the centroids ("coords") of each connected 
%component or region (elements of "CC") in the image and reconfigures the
%structure to an Nx3 array ("centroids").
coords = regionprops(CC, 'Centroid');
centroids = cat(1, coords.Centroid);


end

