%% System Properties
camPixelPitch = 3.75e-6;  % [m]
lensDiameter = 2.54e-2; % [m]
SLMscreenIndex = 3;
LyotCameraID = 2;

%% Lyot Plane Camera
camera = 'DMK42BUC03';
exposure = 1/500;
fps = []; % Not possible to modify for this camera
videoFormat = 'Y800 (1280x960)';
[vid,src] = f_selectCamera(LyotCameraID,camera,exposure,fps,videoFormat);

%% Find center coordinates
[maskShiftY,maskShiftX,systemPupilPixelSIze,mainDataCenter, ...
mainDataRadius] = f_selfCenterSLMmask(camPixelPitch,lensDiameter, ...
SLMscreenIndex,vid,'data');
%% Test shifting
TC = 2; enablechange = false;
[Xr,Yr,aspectRatio,monitorSize] = f_MakeScreenCoords(SLMscreenIndex,false);
% From relative (origin in the center of the screen) to absolute (origin in the upper left corner of the screen):
[maskShiftXsgn, maskShiftYsgn] = f_calcHScoorToSgnCoor(maskShiftX/monitorSize(1), maskShiftY/(monitorSize(2)));
[angularTranstion, rhoRadius] = cart2pol(Xr - maskShiftXsgn, (Yr - maskShiftYsgn)/aspectRatio);
maskGen = @(rho,theta,r0,TC) mat2gray(rho <= r0).*angle(exp(1i*TC*(theta + pi/TC)));
maskr = maskGen(rhoRadius,angularTranstion,0.25/aspectRatio,TC);
imagesc(maskr);