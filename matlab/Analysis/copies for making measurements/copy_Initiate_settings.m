% % % % Initiate camera and rotation stage

% load('Instruments')
% 
% return

%% Camera
Get_camera_on % Tries maximum 500 times to detect
% %  display connected cameras
% % imaqhwinfo
% imaqreset
% % imaqhwinfo
% hw = imaqhwinfo
% % create video input source
% vid = videoinput(hw.InstalledAdaptors{1}, 1, 'Mono14');
% % vid = videoinput('gentl', 1, 'Mono14');
% % src = getselectedsource(vid);

% A = 0;
% return

%% get some hardware info
% get(vid)
% vid_info = imaqhwinfo(vid)
% % vid_info = imaqhwinfo('gentl')
% vid_info = imaqhwinfo('gentl',1)


%% Rotation stage port definition

% %%create serial connection or restore it
% Find a serial port object.
if windows_or_linux == 1 % Windows(0) or Linux(1)
    s = instrfind('Type', 'serial', 'Port', '/dev/ttyUSB0', 'Tag', '');
else
    s = serial('COM5'); % create serial at COMx
end

% Create the serial port object if it does not exist
% otherwise use the object that was found.
if isempty(s)
    s = serial('/dev/ttyUSB0');
else
    fclose(s);
    s = s(1);
end

% return

%% Serial communication with the rotation stage
% % create serial communication channel with EPS300 controler
% serialCom=serial(s);
% communication channel parameters
set(s,'BaudRate',57600);       
set(s,'DataBits',8);            
set(s,'Parity','none');          
set(s,'StopBits',1);         
set(s,'FlowControl','software');     
set(s,'Terminator','CR/LF');
% %open the serial communication channel
% fopen(serialCom);


