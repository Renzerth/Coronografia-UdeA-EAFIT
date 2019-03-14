function [DatalogDir,ProcessedDir,ltcvect,lglvect] = ...
f_CreateFoldersRegisters(maskName,tcvect,glvect,slm,dataDir,outDir,pathSep,infoDelim,dirDelimiter,meas,proc)

%% General saving registers
ltcvect = length(tcvect); % Length of the tc vector
lglvect = length(glvect); % Length of the gl vector
strDate = date; % Today's date is retrieved from the local machine
MeasSize = strcat(maskName,infoDelim,'mask',infoDelim,'tcs',infoDelim, ...
           num2str(ltcvect),infoDelim,'gls',infoDelim,num2str(lglvect));
% Datalog with the number of measurements for tc's and gl's

%% Measurement folder creation (Datalog)
if meas
  % Folder name:
  Datalogfldr = strcat(strDate,infoDelim,slm,infoDelim,MeasSize); 
  % Specific measurement folder:
  f_createNextFolderName(dataDir,Datalogfldr,dirDelimiter);       
end
DatalogDir = strcat(dataDir,pathSep,Datalogfldr); 

%% Output folder creation (processed images)
if proc 
  % Folder name:
  Procfldr = strcat(date,infoDelim,'processed',infoDelim,slm, ...
                       infoDelim,MeasSize); 
  % Specific processing folder:
  f_createNextFolderName(outDir,Procfldr,dirDelimiter); 
end
ProcessedDir = strcat(outDir,pathSep,Procfldr); % Specific processing
                                                     % folder
end 