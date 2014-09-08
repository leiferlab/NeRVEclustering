function output = stackLooper(stack_list)
%STACKLOOPER takes a cell array of .tif 3D image stacks and outputs the
%number of neurons in the image, a 2D Nx3 array of the coordinates of the
%centroids of those neurons, and a slice browser window of the 3D, binary
%image of the neurons. 

output = cell(length(stack_list),1);

    parfor i = 1:length(stack_list)      %for each image stack...

        %Convert the .tif image stack to a 3D volume matrix
        stack_matrix = stackLoad(stack_list{i});
        imsize = size(stack_matrix);
        stack_matrix = image_resize(stack_matrix,imsize(1),imsize(2),imsize(3)*3.45);

        %Convert the 3D matrix representation of the image stack to a 3D matrix
        %representation of the binary image stack with parameters set in fields
        %of "options"

        options = struct('minObjSize', 100, 'maxObjSize', 200, 'minSphericity', 0.9, 'watershedFilter', 1, 'filterSize', [30,30,15], 'noise', 0.5);

%         options.minObjSize = 100;
%         options.maxObjSize = 200;
%         options.minSphericity = 0.9;
%         options.watershedFilter = 1;
%         options.filterSize = [30,30,15];
%         options.noise = .5;
        neurons = WormSegmentHessian3d(stack_matrix, options);

        %Output the coordinates of the centroids ("centroids") of the neurons
        %in the binary image stack ("neurons") and the number of neurons
        %("numNeurons") in the image stack. 
        [centroids, numNeurons] = findCentroids(neurons);
        output{i,1} = centroids;
        %output{i,2} = numNeurons;

        %Open an instance of slicebrowser to convert the 3D matrix
        %representation of the binary image to a viewable image stack.
        %coords = SliceBrowser(neurons);
        i

    end


end