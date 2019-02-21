%% General saving registers
ltcvect = length(tcvect); % Length of the tc vector
lglvect = length(glvect); % Length of the gl vector
strDate = date; % Today's date is retrieved from the local machine
MeasSize = strcat(maskName,'_mask_tcs_',num2str(ltcvect),'_gls_',num2str(lglvect));
% Datalog with the number of measurements for tc's and gl's

%% Measurement folder creation (Datalog)
Datalogfldr = strcat(date,'_',slm,'_',MeasSize); % Folder name
DatalogDir = strcat(dataDir,pathSep,Datalogfldr); % Specific measurement
                                                  % folder
ax = exist(DatalogDir, 'dir'); % 7 if folder exists, 0 if not
if ax ~= 7 % Create a folder if it doesn't exist
    mkdir(DatalogDir);
end

%% Output folder creation (processed images)
Outfldrname = strcat(date,'_processed_',slm,'_',MeasSize); % Folder name
ProcessedDir = strcat(outDir,pathSep,Outfldrname); % Specific processing
                                                   % folder
ax = exist(ProcessedDir, 'dir'); % 7 if folder exists, 0 if not
if ax ~= 7 % Create a folder if it doesn't exist
    mkdir(ProcessedDir);
end 