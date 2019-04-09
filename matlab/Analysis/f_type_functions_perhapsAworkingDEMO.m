%% System Properties
camPixelPitch = 3.75e-6;  % [m]
lensDiameter = 2.54e-2; % [m]
SLMscreenIndex = 3;
LyotCameraID = 2;

%% Lyot Plane Camera
camera = 'DMK42BUC03';
exposure = 1/225;
fps = []; % Not possible to modify for this camera
videoFormat = 'Y800 (1280x960)';
[vid,src] = f_selectCamera(LyotCameraID,camera,exposure,fps,videoFormat);

%% Find center coordinates
[maskShiftY,maskShiftX,systemPupilPixelSIze,mainDataCenter, ...
mainDataRadius] = f_selfCenterSLMmask(camPixelPitch,lensDiameter, ...
SLMscreenIndex,vid,'data');
[maskShiftX, maskShiftY] = f_calcHScoorToSgnCoor(maskShiftX/monitorSize(1), maskShiftY/monitorSize(2));
%% Test shifting
TC = 2; enablechange = false;
[Xr,Yr,aspectRatio,monitorSize] = f_MakeScreenCoords(3,false);
scaledY = Yr/aspectRatio;
[angularTranstion,rhoRadius] = cart2pol(Xr - maskShiftX, scaledY - maskShiftY/aspectRatio);
maskr = mat2gray((rhoRadius <= 0.25/aspectRatio).* ...
                 angle(exp(1i*TC*(angularTranstion + pi/TC))));
imagesc(maskr);