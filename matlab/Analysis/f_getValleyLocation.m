function [rowCoord,colCoord] = f_getValleyLocation(dataMatrix,mode)
%% Find Valley Referece points
cumulativeAverageX = mean(dataMatrix,2);
cumulativeAverageY = mean(dataMatrix,1);
[~,peakLocationX] = findpeaks(cumulativeAverageX,'MinPeakProminence',0.5,'MinPeakHeight',0.75*max(cumulativeAverageX),'NPeaks',2,'MinPeakWidth',20);
[~,peakLocationY] = findpeaks(cumulativeAverageY,'MinPeakProminence',0.5,'MinPeakHeight',0.75*max(cumulativeAverageY),'NPeaks',2,'MinPeakWidth',20);
%% Mode selection
switch (lower(mode))
  case 'vortex'
    if length(peakLocationX) ~= 2 && length(peakLocationY) ~= 2
      warning('Failed to locate peaks.');
      figure('units','normalized','outerposition',[0 0 1 1],'color','white');
      imagesc(dataMatrix);
      rawCoordinates = round(ginput(1));
      rowCoord = rawCoordinates(2);
      colCoord = rawCoordinates(1);
    else
      peakLocationX = sort(peakLocationX,'ascend');
      peakLocationY = sort(peakLocationY,'ascend');
      
      valleyPointX = min(cumulativeAverageX(peakLocationX(1):peakLocationX(2)));
      valleyPointY = min(cumulativeAverageY(peakLocationY(1):peakLocationY(2)));
      
      [rowCoord,~] = ind2sub(size(cumulativeAverageX),find(cumulativeAverageX == valleyPointX)); % Find valley point coordinates
      [~,colCoord] = ind2sub(size(cumulativeAverageY),find(cumulativeAverageY == valleyPointY));
    end
  case 'single'
    if length(peakLocationX) ~= 2
      peakLocationY = sort(peakLocationY,'ascend');
      valleyPointY = min(cumulativeAverageY(peakLocationY(1):peakLocationY(2)));
      [~,colCoord] = ind2sub(size(cumulativeAverageY),find(cumulativeAverageY == valleyPointY));
      rowCoord = {};
    elseif length(peakLocationY) ~= 2
      peakLocationX = sort(peakLocationX,'ascend');
      valleyPointX = min(cumulativeAverageX(peakLocationX(1):peakLocationX(2)));
      [rowCoord,~] = ind2sub(size(cumulativeAverageX),find(cumulativeAverageX == valleyPointX));
      colCoord = {};
    end
  otherwise
    error('Invalid Option.');
end
end