%% General saving registers
ltcvect = length(tcvect); % Length of the tc vector
lglvect = length(glvect); % Length of the gl vector
strDate = date; % Today's date is retrieved from the local machine
MeasSize = [maskName '_mask_tcs_' num2str(ltcvect) '_gls_' num2str(lglvect)];
% Datalog with the number of measurements for tc's and gl's
switch windows
    case 7
        pathSep = '/'; % Path separator  
    case 10
        pathSep = '\'; % Path separator  
end

%% Measurement folder creation (Datalog)
Datalogfldr = [date '_' MeasSize]; % Folder name
cd(dataDir);
ax = exist(Datalogfldr, 'dir'); % 7 if folder exists, 0 if not
if ax ~= 7 % Create a folder if it doesn't exist
    mkdir(Datalogfldr);
end
DatalogDir = [dataDir pathSep Datalogfldr]; % Specific measurement folder

%% Output folder creation (processed images)
Outfldr = [date '_processed_' MeasSize]; % Folder name
cd(outDir); % Goes to the output directory
ax = exist(Outfldr, 'dir'); % 7 if folder exists, 0 if not
if ax ~= 7 % Create a folder if it doesn't exist
    mkdir(Outfldr);
end
ProcessedDir = [outDir pathSep Outfldr]; % Specific measurement folder
cd(analysDir);