function [dataSize,relativeCoord] = f_computeCenterFromImageMin()
%% Load data
[centerLinesData,~,~,~] = f_Mod_file_loader('.jpg',0);
%% Compute Center Point
[rowCoord,~] = f_getValleyLocation(centerLinesData{1},'single');
[~,colCoord] = f_getValleyLocation(centerLinesData{2},'single');
centerPoint = [rowCoord(1),colCoord(2)];
%% relative coordinates
dataSize = size(centerLinesData{1});
relativeCoord = centerPoint./dataSize;
end