function [imgpartPath,measfullpath,ProcessedDir,ltcvect,lglvect,totalImgs,numberedFolderMeas] = ...
f_CreateFoldersRegisters(maskName,tcvect,glvect,slm,cameraPlane,dataDir,...
outDir,pathSep,infoDelim,dirDelim,lastmeas,meas,proc)

%% General saving registers
ltcvect = length(tcvect); % Length of the tc vector
lglvect = length(glvect); % Length of the gl vector
totalImgs = ltcvect*lglvect; % Number of images to be taken
strDate =date; % Today's date is retrieved from the local machine
MeasSize = strcat(maskName,infoDelim,'mask',infoDelim,'tcs',infoDelim, ...
           num2str(ltcvect),infoDelim,'gls',infoDelim,num2str(lglvect));
% Datalog with the number of measurements for tc's and gl's

%% Measurement folder creation (Datalog)
%%% Last measurement's date directory
lastmeasDir = strcat(dataDir,pathSep,lastmeas);

%%% Asks if a measurement will be performed
if meas
  createFoldMeas = 1;
  %%% Folder name:
  Measfldr = strcat(strDate,infoDelim,slm,infoDelim,MeasSize);  % Use the wanted measurement's name
  % Explanation: save(directory,filename,variables) % ,'-append'
  save(lastmeasDir,'Measfldr'); % Save as .mat
  
else % meas == 0
  createFoldMeas = 0;
  struct = load(strcat(lastmeasDir,'.mat'));  
  Measfldr = struct.Measfldr; % Use last measurement's name
  %%% Folder name:
   
end
%%% Specific measurement folder:
numberedFolderMeas = f_createNextFolderName(dataDir,Measfldr, ...
                                          dirDelim,pathSep,createFoldMeas);                        
MeasDir = strcat(dataDir,pathSep,numberedFolderMeas); % Datalog folder that
% takes into account if the folder alraedy exists (numbered folder)
imgpartPath = strcat(MeasDir,pathSep,cameraPlane,infoDelim); % Partial path. More 
%%% Path of the cell with all the measurements
measfullpath = strcat(imgpartPath,'allmeas'); % Name of the saved cell of 
% images. More information will be concatenated for a full path of the 
% measured images inside the two "for" loops in f_AutomateMeasurement. As 
% well, the variable "imgfullpath" is created inside this function

%% Output folder creation (processed images)
% Folder name:
Procfldr = strcat(strDate,infoDelim,'processed',infoDelim,slm, ...
                  infoDelim,MeasSize); 
if proc 
  createFoldProc = 1;
else
  createFoldProc = 0;
end
% Specific processing folder:
numberedFolderProc = f_createNextFolderName(outDir,Procfldr,dirDelim, ...
                                            pathSep,createFoldProc); 
ProcessedDir = strcat(outDir,pathSep,numberedFolderProc); 
end 