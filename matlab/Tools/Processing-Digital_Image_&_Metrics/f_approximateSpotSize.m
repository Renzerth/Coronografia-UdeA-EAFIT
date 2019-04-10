function [x,y,regionCentroid,aproxRadius] = f_approximateSpotSize(cdata,varargin)
% Approximates a PSF radius by binarizing the image and finding its radius by two methods:
% blobRadius and radiusArea -> averages these results
%
% The input image could be circular or somehow different, for a function that demands the 
% shape to be circular, use "f_findCircleShapedIntensity.m"
%
% Inputs:
% cdta: input image (2D matrix, already with the " im2double" previously applied)
% threshold: optional number in [0,1] for the binarization. It is a value
% of sensitivity and is not directly associated to the intensity of the
% image
%
% Outputs:
% x,y: vectors with the lamda/D scaling
% regionCentroid: coordinates of the spot in the [row, col] = [y,x] format
% and in approximated pixels (integer)
% aproxRadius: approximate spot radius

%% Threshold as input or default
if nargin == 2
    Threshold = varargin{1}; % User input
else
    Threshold = 0.3; % 70 percent of the energy  is considered. One only analyzes this energy since it is significative in
                               % an Airy pattern. It is assumed that the central intensity defines the spot's axis or 
end
[sizeY,sizeX] = size(cdata);
binaryData = false(sizeY,sizeX);
binaryData(cdata>=Threshold) = true; % (1-Threshold) percent of the energy
% OLD binaryData = im2bw(cdata,Threshold); % use imbinarize for newer MATLAB versions (2017 and further versions)

%% Estimation of the pixel size of the Airy Disk 
regionInfo = regionprops(binaryData,'Centroid','area','MajorAxisLength','MinorAxisLength');
[~, sortIndexes] = sort(cat(1,regionInfo.Area), 'descend');
mainIndex = sortIndexes(1);

%% Centroid of the spot
regionCentroid = round(fliplr(cat(1,regionInfo(mainIndex).Centroid))); % Row, Column format coordinates in approximated pixels

%% Radius estimative between the major and the minor axis of the elliptical spot
blobRadius = mean([regionInfo(mainIndex).MajorAxisLength, regionInfo(mainIndex).MinorAxisLength],2)/2;

%% Morphological radius estimative from the total area of the spot: it assumes it is a perfect circle
radiusArea = sqrt(max(regionInfo(mainIndex).Area)/pi); 
aproxRadius = ceil(mean([blobRadius,radiusArea])); % ceil once a pixel is occupied, it is assumed to fill the whole pixel

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