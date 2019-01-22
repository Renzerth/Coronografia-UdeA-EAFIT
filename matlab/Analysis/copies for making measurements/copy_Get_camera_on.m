

try
imaqreset
hw = imaqhwinfo;
vid = videoinput('gentl', 1, 'Mono14');  % create video input source
% hw.InstalledAdaptors{1}
vid_device = getselectedsource(vid);

vid_device.Gain = 0; % Gain is digital and induces non-linear effects
vid_device.ExposureTime = 1000; %50000;ge
vid.TriggerRepeat = 0;
vid.FramesPerTrigger = 4; % Default frames to capture per trigger

preview(vid)

catch
  Get_camera_on
end


