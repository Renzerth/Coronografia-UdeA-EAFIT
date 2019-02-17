%% Initialization
experimentValues = [05,0,0]; % [PH,TC,DST] Used during measurements
totalImages = 75;
experimentalImages = cell(1,totalImages);
index = 1;
previousIndex = 1;
errorCount = 0;
%% Set experiment voltages
minVoltage = 0.6;
maxVoltage = 75; % This voltage value cannot suprass device allowed max voltage
voltageStep = (maxVoltage-minVoltage)/totalImages;
voltageRange = minVoltage:voltageStep:maxVoltage-voltageStep;
%% Set Camera
camera = 'DMK42BUC03';
exposure = 1.282;
videoFormat = 'Y800 (1280x960)';
[vid,src] = CAMERA.selectCamera(camera,exposure,videoFormat);
recordingDelay = 1.5;
%% Set serial channel
COMPort = 'COM4';
[serialObject] = TOOLS.setPiezoSerialChannel(COMPort);
terminator = '\r\n';
fopen(serialObject);
%% Set save directory
if experimentValues(2) ~= 0
    targetFolder = sprintf('ILPC_VOTC%i_Interf_PH%i',experimentValues(2),experimentValues(1));
else
    targetFolder = sprintf('ILPC_REF_Interf_PH%i',experimentValues(1));
end
[folderPath] = TOOLS.makeParentFolder(0,targetFolder);
fileFormat = '.png';
%% Image Tag
[imageName] = TOOLS.makeImageTag('PCSVortexInterferograms');
%% Capture and displacement
while index < totalImages
    try
        for index = previousIndex:totalImages
            fprintf(serialObject,sprintf('XV%1.1f\n',voltageRange(index))); % Send voltage value to piezo driver
            pause(recordingDelay);
            SingleFrame = getsnapshot(vid);
            experimentalImages{index} = SingleFrame;
            fileName = sprintf(imageName,experimentValues,voltageRange(index));
            imwrite(SingleFrame,strcat(folderPath,'\',fileName,fileFormat));
        end
    catch ME
        index = previousIndex;
        errorCount = errorCount + 1;
        if errorCount > 5
            error('Failed to accomplish the task due to unavailable COM ports or locked Camera.');
        end
        fprintf('Unable to proceed. -- Attempt:%i -- Reason:\n',errorCount);
        fprintf(strcat(ME.message,'\n'));
        fprintf('Close all other camera programs or verify COM port.\n');
        fprintf('System will try again with procedure.\n');
        if strcmp(serialObject.Status,'closed')
            fopen(serialObject);
        end
        pause(3);
    end
end
%% Set to zero value
fprintf(serialObject,'XV0.0\n');
%% End notification
for beepTimes = 1:3
    beep();
    pause(0.2);
end
%% Close port
fclose(serialObject);
delete(vid);