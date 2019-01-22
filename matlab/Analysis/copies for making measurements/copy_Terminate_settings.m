% Terminate_settings

%% Rotation stage

s = instrfind;
disp(s);
fclose(s);
fopen(s);
% resetPos(s,0)
fclose(s);
% %find open communication devices, close them
% openDevices=instrfind('Status','open');
% for i=1:length(openDevices)
%     fclose(openDevices(i));   
% end

delete(s)
clear s
instrfind

%% Camera
delete(vid)
clear all, close all;


