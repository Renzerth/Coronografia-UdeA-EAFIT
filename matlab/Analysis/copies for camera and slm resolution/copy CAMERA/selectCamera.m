function [vid,src] = selectCamera(camera,exposure,format)
%% Select camera for predefined settings
try
switch camera
  case 'DMK42BUC03'
    vid = videoinput('tisimaq_r2013', 1, format);
    src = getselectedsource(vid);
    src.Brightness = 0;
    src.Contrast = 0;
    src.Exposure = exposure; %PH50um : 0.0030 -- Y0.0256 :: PH25um : 0.0083 -- Y0.0556  :: PH15um --Y0.1429 : : 0.0227 :: PH10um : 0.0435 -- Y0.200 PH5um : 0.363 -- Y3.099
    src.GainAuto = 'Off';
    src.Gain = 34;
    
  case 'DMK41BU02.H'
    vid = videoinput('tisimaq_r2013', 1, format);
    src = getselectedsource(vid);
    src.Exposure = exposure;
    src.GainAuto = 'Off';
    
  case 'DMK23U445'
    vid = videoinput('tisimaq_r2013', 1, format);
    src = getselectedsource(vid);
    src.Exposure = exposure;
    src.GainAuto = 'Off';
    
  otherwise
    warning('Undefined Camera');
end
catch ME
    error('Camera is connected but properties cannot be accessed. Verify the USB Port.');
end