function [systemPupilPixelSIze] = computePupilPixelSize(spotPixelSize, detectorPixelPitch,systemPupilSIze)
% systemPupilSIze, apRad IN CM

%% Linear relationship
spotDiameter = 2*spotPixelSize;
physicalSize = spotDiameter*detectorPixelPitch;
illuminationRatio = physicalSize/systemPupilSIze;
systemPupilPixelSIze = round(spotDiameter/illuminationRatio);
end