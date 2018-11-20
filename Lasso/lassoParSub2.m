function  lassoParSub2(dataFolder)
%%set up data
numTrial = 10;

fileName = [dataFolder filesep 'fixed_data.mat'];
load(fileName)
newRatio2 = newNA;
mode3 = smooth(behaviorZ,20);
dtheta = eigAngleFilt;
nNeuron = length(newRatio2(:,1));
if nNeuron==0
    smoothNA = 0;
else 
    for i =1:nNeuron
        smoothNA(i,:) = smooth(newRatio2(i,:),10); 
    end
end

cutPoint = floor(.6*length(eigAngleFilt));
%%

binNum = 5;
goalPts = ceil((.1*cutPoint)/binNum);
train = eigAngleFilt(1:cutPoint);

h_train = histogram(train,binNum);


edgesTrain = h_train.BinEdges;

[minCount, minIdx] = min(h_train.Values);

if minCount < goalPts
    binNum = 4;
    h_train = histogram(train,binNum);
end
[minCount, minIdx] = min(h_train.Values);
edgesTrain = h_train.BinEdges;

percPick = goalPts/minCount;
numPts = ceil(percPick*minCount);

for i = 1:binNum
   bin(i).allPtsTrain = find(train>edgesTrain(i) & train< edgesTrain(i+1));
end
%%
for j = 1:numTrial

for i = 1:binNum
    temp = randperm(length(bin(i).allPtsTrain));
   binPts(i,:) = bin(i).allPtsTrain(temp(1:numPts));
end

trainPts(j,:)  = sort(reshape(binPts,[1, binNum*numPts]));

end
%%

parfor k = 1:numTrial

training = trainPts(k,:);


alpha = 0.001:.5:1;
alpha = [alpha 1];
timeStep = 1:10;

for j = 1:length(alpha)
    for dt = 1:length(timeStep)
        [fold(k).alph(j).regress(dt).WeightsdTheta,fold(k).alph(j).regress(dt).FitdTheta,fold(k).alph(j).regress(dt).WeightsMode3,fold(k).alph(j).regress(dt).FitMode3, fold(k).alph(j).regress(dt).DataMat]   = lassoPickPts(smoothNA, timeStep(dt), dtheta, mode3,  alpha(j),training);
    
    end
end

end
save([dataFolder filesep 'lassoOutput2.mat'], 'fold', '-v7.3');