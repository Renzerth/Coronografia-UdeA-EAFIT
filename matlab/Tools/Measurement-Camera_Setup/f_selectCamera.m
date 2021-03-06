function [vid,src,fps] = f_selectCamera(cameraID,camera, ...
                                              exposure,fps,format)
% Inputs:
%   camera
%   exposure
%   fps
%   format
%
% Outputs:
%   vid: video input object
%   src: currently selected video source object

%% Select camera for predefined settings
try
switch camera
  case 'DMK42BUC03' % CMOS camera used for the transmission SLM
    vid = videoinput('tisimaq_r2013', cameraID, format);
    src = getselectedsource(vid);
    src.Brightness = 0;
    src.Contrast = 0;
    src.Exposure = exposure; %PH50um : 0.0030 -- Y0.0256 :: PH25um : 0.0083 -- Y0.0556  :: PH15um --Y0.1429 : : 0.0227 :: PH10um : 0.0435 -- Y0.200 PH5um : 0.363 -- Y3.099
%     src.FrameRate = fps; Not Available
    fps = '25'; % fps: not possible to change for this camera
    src.GainAuto = 'Off'; % Gain auto feature is always off
%   src.Gain = 0; % In: [0,30]
    
 case 'DMK23U445' % PSF plane % CMOS
    vid = videoinput('tisimaq_r2013', cameraID, format);
    src = getselectedsource(vid);
    src.Exposure = exposure;
    src.FrameRate = fps; % set the fps
    fps = src.FrameRate; % read the fps
    src.GainAuto = 'Off'; % Gain auto feature is always off

  case 'DMK41BU02.H'
    vid = videoinput('tisimaq_r2013', cameraID, format);
    src = getselectedsource(vid);
    src.Exposure = exposure;
    src.FrameRate = fps; % set the fps
    fps = src.FrameRate; % read the fps
    src.GainAuto = 'Off'; % Gain auto feature is always off
    
  otherwise
    warning('Undefined Camera');
end
catch ME
    error('Camera is connected but properties cannot be accessed. Verify the USB Port.');
end

end