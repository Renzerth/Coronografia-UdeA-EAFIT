function [systemPupilPixelSize] = f_computePupilPixelSize(spotPixelSize, detectorPixelPitch,systemPupilSize)

% Inputs:
%   spotPixelSize: radius of the spot in pixels (detector plane)
%   detectorPixelPitch:  in [m/pix]
%   systemPupilSIze: in [m]
%
% Outputs:
%   systemPupilPixelSIze in pixels

%% Linear relationship between the spotPixel size and the system's pupil size
spotDiameter = 2*spotPixelSize; % Double of the radius
physicalSize = spotDiameter*detectorPixelPitch; % in meters
illuminationRatio = physicalSize/systemPupilSize; % A proportion
systemPupilPixelSize = round(spotDiameter/illuminationRatio); % nearest value in pixels
end