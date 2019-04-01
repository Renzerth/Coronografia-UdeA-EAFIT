function [x,y,regionCentroid,aproxRadius] = f_approximateSpotSize(cdata,varargin)
% Approximates a PSF radius by binarizing the image and finding its radius by two methods:
% blobRadius and radiusArea -> averages these results
%
% Inputs:
% cdta: input image (2D matrix, already with the " im2double" previously applied)
%
%
% Outputs:
% x,y: vectors with the lamda/D scaling
% regionCentroid: coordinates of the spot
% aproxRadius: approximate spot radius

%% Threshold as input or default
if nargin == 2
    Threshold = varargin{1}; % User input
else
    Threshold = 0.01; % From 0 to 1
end

%% Estimation of the pixel size of the Airy Disk 
binaryData = im2bw(cdata,Threshold); % use imbinarize for newer MATLAB versions (2017 and further versions)
regionInfo = regionprops(binaryData,'Centroid','area','MajorAxisLength','MinorAxisLength');
[~, sortIndexes] = sort(cat(1,regionInfo.Area), 'descend');
mainIndex = sortIndexes(1);

%% Centroid of the spot
regionCentroid = cat(1,regionInfo(mainIndex).Centroid);

%% Radius stimative between the major and the minor axis of the elliptical spot
blobRadius = mean([regionInfo(mainIndex).MajorAxisLength, regionInfo(mainIndex).MinorAxisLength],2)/2;

%% Morphological radius estimative from the total area of the spot: it assumes it is a perfect circle
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
figure;
imagesc(x,y,cdata);
title('Test for f-approximateSpotSize.m')
end