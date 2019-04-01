function [x,y,regionCentroid,aproxRadius] = f_approximateSpotSize(cdata,varargin)

%% Threshold as input or default
if nargin == 2
    Threshold = varargin{1}; % User input
else
    Threshold = 0.01; % From 0 to 1
end

%% Estimation of the pixel size of the Airy Disk 
binaryData = im2bw(cdata,Threshold); % use imbinarize for newer MATLAB versions (2018 and further versions)
regionInfo = regionprops(binaryData,'Centroid','area','MajorAxisLength','MinorAxisLength');
[~, sortIndexes] = sort(cat(1,regionInfo.Area), 'descend');
mainIndex = sortIndexes(1);
blobRadius = mean([regionInfo(mainIndex).MajorAxisLength, regionInfo(mainIndex).MinorAxisLength],2)/2;
regionCentroid = cat(1,regionInfo(mainIndex).Centroid);
radiusArea = sqrt(max(regionInfo(mainIndex).Area)/pi); 
aproxRadius = ceil(mean([blobRadius,radiusArea])); % Once a pixel is occupied, it is assumed to fill the whole pixel

%% Airy Disk - Image size ratio
[p,q] = size(cdata);
halfP = p/2;
halfQ = q/2;
xLengthRatio = halfP/aproxRadius;
yLengthRatio = halfQ/aproxRadius;

%% Scalled coordinates
x = -xLengthRatio:xLengthRatio-1;
y = -yLengthRatio:yLengthRatio-1;

%% Plot
imagesc(x,y,cdata);
end