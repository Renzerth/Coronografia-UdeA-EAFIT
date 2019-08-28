function [mainDataCenter,mainDataRadius,dataProportion] = ...
                              f_findCircleShapedIntensity(rawData,varargin)
% The input image is assumed to be circular, for a more general function,
% use "f_approximateSpotSize.m"

 % Inputs:
%   rawData: double format 2D image (single matrix)
%   plotsOn: true or false
%
% Outputs:
%   mainDataCenter: 
%   mainDataRadius: in pixels
%   dataProportion:

%% Settings
if nargin == 2
    plotsOn = varargin{1};
else
    plotsOn = true;
end

%% Parameters
luminanceThresh = 0.1; %0.1
radiusRangePerce = 0.15; %0.25

%% Process Reference
shapeData = rawData-mean(rawData(:)); % Background filtering
enhancedRange = mat2gray(log10(im2double(shapeData + 1))); % avoid -Inf log
                                                % and turn it back to uint8

%% Detect data regions
totalCounts = numel(enhancedRange);
% binaryData = imbinarize(enhancedRange,luminanceThresh);
binaryData = double(enhancedRange>luminanceThresh);
binaryData = imfill(binaryData,4,'holes');
binaryData = imopen(binaryData,strel('disk',6)); % Binary shape must fill
                         % circle to increase enclosed area approximation

%% Compute blob properties for the pixel size of the circle
regionInfo = regionprops(binaryData,'Centroid','area', ...
                                    'MajorAxisLength','MinorAxisLength');
[~, sortIndexes] = sort(cat(1,regionInfo.Area), 'descend');
if isempty(sortIndexes)
  error('No circular intensity was registered by the camera.');
end
mainIndex = sortIndexes(1);
blobRadius = mean([regionInfo(mainIndex).MajorAxisLength, ...
                   regionInfo(mainIndex).MinorAxisLength],2)/2;
regionCentroid = cat(1,regionInfo(mainIndex).Centroid);

%% Compute area approximation properties
areaCounts = histcounts(binaryData);
dataProportion = areaCounts/totalCounts;
radiusEstimation = sqrt(areaCounts(2)/pi);
dataRanges = round([1-radiusRangePerce,1+radiusRangePerce]* ...
                   mean([radiusEstimation,blobRadius]));

%% Generate circle with the Circle Hough Transform (CHT)
% dataRanges (the mean radius) is an input for using the CHT.
[dataCenter, dataRadii, circMetrics] = imfindcircles(binaryData, ...
                                              dataRanges,'Sensitivity', 1);
[~,valueIndex] = max(circMetrics);

%% Mean of the centroid (blob) and the center (CHT) methods
mainDataCenter = mean([dataCenter(valueIndex,:); regionCentroid],1);

%% Mean of the blob radius (region major/minor axis estimation) and the CHT radius (area of a circle assumed)
mainDataRadius = mean([dataRadii(valueIndex,:),blobRadius]);

%% Validate wether one has a circular region
if abs(diff([dataRadii(valueIndex,:),blobRadius]))/blobRadius > 0.06
    error('Blob region is not a circular area'); % Compare the difference
                         % from circular prediction to the blob behaviour
end

%% Visualize detected data
if plotsOn
    imagesc(rawData);
    viscircles(mainDataCenter,mainDataRadius,'color','b');
end

end