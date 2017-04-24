function output=distanceClassify(centers,samples,hitCutoff)
%%%
% distanceClassify matches each row of samples to the row in centers which
% it is closest to. Each center can only be represented by at most one
% sample. The program works by calculating the correlation between
% each row of sample and each row of centers, and then matches based on the
% largest correlation. Matches are only assigned of the correlation is
% larger that hitCutoff

%INPUTS:
%   centers - c x n matrix where each row is the center of a
%   classification group in n dimensinos
%
%   samples - k x n matrix where each row is a sample point to be matched
%   with one of the centers.
%
%   hitCutoff - the minimum correlation needed before samples are matched
%   with a center

%OUTPUT:
%   output - a cell with
%       output{1} - k x 1 vector with the center each row of samples is
%       matched too. Unmatched rows are given an Nan
%       output{2} - k x 1 vector with the correlation between each row of
%       samples and the center it was matched to. Unmatched rows are given
%       0. 

if nargin==2
hitCutoff=.2;
end

if isempty(samples)
    id=[];
    weights=[];
    output={id,weights};
    return
elseif size(samples,2)~= size(centers,2)
    error(' Dimension of centers and samples must match!')
end
samples=bsxfun(@rdivide,samples,sqrt(sum(samples.^2,2)));
matchProjectionsRaw=centers*samples';
matchProjections=matchProjectionsRaw.*(matchProjectionsRaw>hitCutoff);

aMax = max(matchProjections);
%if there is no hit, set the max high, the projection will never equal it
aMax(aMax==0)=10;
%find where the best match occurs
matchBool=bsxfun(@eq, matchProjections,aMax);
matchBool=matchProjections.*matchBool;
[match_x,match_y,val]=find(matchBool);


%% clear doubles
%find if there are time points with doubles
doubles=check4doubles(match_x);


    for iiCheck=doubles'
        % get the neurons that were assigned to the same index
        doubleSearch=find(match_x==iiCheck);
        % get the weights for each of them
        w=val(doubleSearch);
        % remove the one with the smaller projection. 
        idx2Remove=doubleSearch(w~=max(w));
        match_x(idx2Remove)=nan;
        % also set that weight to zero
        val(idx2Remove)=0;
    end

id= nan(size(samples,1),1);
id(match_y)=match_x;
weight=zeros(size(samples,1),1);
weight(match_y)=val;
output={id,weight};

