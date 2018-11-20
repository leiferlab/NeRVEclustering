function [lassoWb1, lassoFitb1,lassoWb2, lassoFitb2,  timeLagMatrix] = lassoPickPts(smoothNA,dt, behav1, behav2, alpha, training)
nNeuron = length(smoothNA(:,1));
dataMat = smoothNA(:,training+99);
timeLagMatrix = smoothNA(:,100-dt:end-100);
for i =2:dt
   dataMat(length(dataMat(:,1))+1:nNeuron*i,:) =  smoothNA(:,training+99-i);
   startN = length(timeLagMatrix(:,1)) +1;
   timeLagMatrix(startN:startN+nNeuron-1, :) = smoothNA(:,100-dt-i:end-i-100);
end
x = linspace(0,7,100);
lam = exp(x);
lam = lam/400;
[lassoWb1, lassoFitb1] = lasso(dataMat', behav1(training), 'LambdaRatio', 0, 'Alpha', alpha, 'Lambda',lam );
[lassoWb2, lassoFitb2] = lasso(dataMat', behav2(training), 'LambdaRatio', 0,'Alpha', alpha,'Lambda',lam);