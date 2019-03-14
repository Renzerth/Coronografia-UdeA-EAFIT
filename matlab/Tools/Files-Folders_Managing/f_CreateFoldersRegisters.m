function [DatalogDir,numberedFolderMeas,ProcessedDir,numberedFolderProc,ltcvect,lglvect] = ...
f_CreateFoldersRegisters(maskName,tcvect,glvect,slm,dataDir,outDir,pathSep,infoDelim,dirDelim,meas,proc)

%% General saving registers
ltcvect = length(tcvect); % Length of the tc vector
lglvect = length(glvect); % Length of the gl vector
strDate = date; % Today's date is retrieved from the local machine
MeasSize = strcat(maskName,infoDelim,'mask',infoDelim,'tcs',infoDelim, ...
           num2str(ltcvect),infoDelim,'gls',infoDelim,num2str(lglvect));
% Datalog with the number of measurements for tc's and gl's

%% Measurement folder creation (Datalog)
% Folder name:
Datalogfldr = strcat(strDate,infoDelim,slm,infoDelim,MeasSize); 
if meas
  % Specific measurement folder:
  numberedFolderMeas= f_createNextFolderName(dataDir,Datalogfldr,dirDelim,pathSep);       
end
DatalogDir = strcat(dataDir,pathSep,Datalogfldr); 

%% Output folder creation (processed images)
% Folder name:
Procfldr = strcat(date,infoDelim,'processed',infoDelim,slm, ...
                     infoDelim,MeasSize); 
if proc 
  % Specific processing folder:
  numberedFolderProc = f_createNextFolderName(outDir,Procfldr,dirDelim,pathSep); 
end
ProcessedDir = strcat(outDir,pathSep,Procfldr); 
end 