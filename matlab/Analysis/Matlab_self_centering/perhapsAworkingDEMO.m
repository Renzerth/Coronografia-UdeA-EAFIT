%% Lyot Plane Camera
camera = 'DMK42BUC03';
exposure = 1/225;
videoFormat = 'Y800 (1280x960)';
[vid,src] = f_selectCamera(2,camera,exposure,videoFormat);
%% System Properties
camPixelPitch = 3.75e-6;  % [m]
lensDiameter = 2.54e-2; % [m]
SLMscreenIndex = 3;
%% Self Centering
[maskShiftY, maskShiftX, systemPupilPixelSIze, mainDataCenter, mainDataRadius] = f_selfCenterSLMmask(camPixelPitch, lensDiameter, 3, vid);