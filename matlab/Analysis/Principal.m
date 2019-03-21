%%  Phase mask SLMs display for the digital vortex coronagraphy and 
%%% free-space simulation of two stars
% Inputs:
% See "Paramters.m"
%
% Outputs:
%  Plots on the PC or on the SLMs
%  Camera captures
%  Processed vortex images (or plots)
%
% Notes:
%  -It is recommended to have your screen's wallpaper as a white color
%  -Units: a.u (arbitrary units) and cm for lengths, radians for angles and
%  um for wavelengths
%  -Helicoidal phase masks are inside the Laguerre-Gauss beams category
%  -The variable mask is complex and is wrappped: mask = exp(i*mask)
%  -All wrappped phases are shown on [-pi,pi] 
%  -Always execute the program whenver you are exactly inside its folder
%  -The Zernike polynomials are always used with a square size and they 
%   stop being almost-orthogonal if the pupil gets truncated (non-circular 
%   anymore)
%
% Bugs:
%  -When the ticks of -pi, ..., pi are shown, sometimes the may move
%  slightly from the real value
%  -The algorithm shall only be executed inside the analysis folder
%  - Algorithm not tested yet in Linux but in Windows 7 and 10
%  -Simulations are not fully accurate with the calculations
%  -Zernike polynomials are not fully wrapped and the "a" and "frac" 
%   constants have not been fully understood
%  -If a camera is not recognized or even the screen index are not acting
%  coherently, restart MATLAB
%  -For high resolution monitors, sometimes the figure bar may appear but
%  this presents no problem for the mask projection
%
% Samuel Plazas(PA1/PA2/TG) - Juan Jose Cadavid(Master thesis) - 2018/2019

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% PHASE MASK GENERATION ON THE SLM's 
%% Parameters and initialization
clc; clear; close all; % Initialization
Parameters; % Adds to the algorithm all the needed parameters
% open Parameters; % Keep open always

%% Directories and add functions
% Adds initial functions
addpath(strcat('..',pathSep,toolsFldr,pathSep,filemanag)); 
[analysDir,toolsDir,dataDir,outDir] = ...
f_addDirectories(analysFldr,toolsFldr,dataFlrd,outFlrd,pathSep);
% Adds all the directories to use in the algorithm

 %% Hardware initialization
 % This is here because the camera is needed for the self-centering
 % algorithm inside "f_DefineSpace"
 if measSimulated == 0 && meas == 1 % When a real measurement will be performed
  % InitializeHardware;
  % Turns the camera on and create all the needed variables. Remember to 
  % leave the preview open
  [vid,~] = f_selectCamera(camera,exposure,format);
 else
     vid = [];
 end

%% Spatial definitions
[rSize,x,y,Xslm,Yslm,rSLM,phiSLM,Xpc,Ypc,rPC,phiPC,monitorSize, ...
shiftCart] = f_DefineSpace(vid,sSupport,sSize,shiftCart,shiftMask,pixSize, ...
scrnIdx,circularMask,z_pupil,coordType,MaxMask,maskSel);
% Defines the cartesian/polar coordinates, its sampling interval and 
% discretized angular part for gl
%%% Coordinates selection:
[X,Y,r,phi] = f_SelectCoordinates(Xslm,Yslm,rSLM,phiSLM,Xpc,Ypc,rPC, ...
                                  phiPC,plotMask);
% f_DefineSpace and f_SelectCoordinates were divided so that in 
% f_SelectCoordinates one just selects already-calculated variables
                
%% Phase mask selection and plot on the screen or on the SLM
% Don't plot the mask if one will process or measure
if proc == 1 || meas == 1
    plotMask = 0;
end
[mask,wrapMask,~,maskName] = f_PlotSelectedMask(X,Y,r,phi, ...
phaseValues,tc,s,ph0,p,WsizeRatio,L,f_FR,bcst,period,T0,frkTyp,Aalpha, ...
Angalp,Angbet,z_coeff,z_a,z_pupil,z_disp_wrap,z_plot,normMag, ...
binMask,binv,MaskPupil,rSize,monitorSize,scrnIdx,coordType,abs_ang, ...
MaxMask,plotMask,maskSel);      
% Dependencies:
% f_PlotSelectedMask -> f_SpiralMask (or any other) -> f_ProjectMask -> 
% f_MaskWrapCircDiscret -> f_discretizeMask -> f_wrapToRange 
%                       -> f_ScaleMatrixData       

              


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% MEASUREMENT/PROCESSING BY AN AUTOMATED PARAMETER VARIATION
%% Folders and register creations on Data and Output    
[imgpartPath,measfullpath,ProcessedDir,ltcvect,lglvect,totalImgs, ...
numberedFolderMeas] = f_CreateFoldersRegisters(maskName,tcvect,glvect, ...
slm,cameraPlane,dataDir,outDir,pathSep,infoDelim,dirDelim,meas,proc);

%% Measurement
if meas

 %% Measurement debugging
 % Usefull for aligning the vortex and adjusting exposure parameters
 if measDebug == 1 && measSimulated == 0
  SingleFrame = f_CaptureImage(vid,dataDir,filename,imgformat,pathSep, ...
                               snapsfldr); 
  % Takes a camera shot,shows a figure and saves it   
  figure; imhist(SingleFrame); % Shows a histogram of the snapshot

  %%% Get some hardware/software/tools info:
  % get(vid): displays the general parameters of the camera
  % getselectedsource(vid): an existing variable (src). Similar to get(vid)
  % imaqhwinfo(vid): displays the driver connection with MATLAB
  % imaqhwinfo: all the installed adaptors in the Image Acquisition Toolbox
  % disp(vid): displays acquisition information
  % imaqtool: toolbox for the camera
  % imaqreset: refresh image acquisition hardware by restoring the settings
 else % Real measurement

  %% Automated measurement
  [refmeasfullpath] = f_AutomateMeasurement(vid,Xslm,Yslm,rSLM,phiSLM, ...
  Xpc,Ypc,rPC,phiPC,s,ph0,p,WsizeRatio,L,f_FR,bcst,period,T0,frkTyp, ...
  Aalpha,Angalp,Angbet,z_coeff,z_a,z_pupil,z_disp_wrap,z_plot,normMag, ...
  binMask,binv,MaskPupil,rSize,monitorSize,scrnIdx,coordType,abs_ang, ...
  MaxMask,maskSel,ltcvect,lglvect,totalImgs,wait,imgpartPath, ...
  dataformat,imgformat,measfullpath,infoDelim,cameraPlane,tcvect, ...
  glvect,measSimulated,recordingDelay,whiteblackref); 
  % Performs measurements and stores them
 end

 %% Termination
 % TerminateSettings; % Future script % Clears variables, closes all and
                       % deactivates the cameras
 % delete(vid); % Clean up the camera

 %% End notification
 if beepSound == 1
     for beepTimes = 1:3 % Numbe of beeps
         beep();
         pause(0.2); % Time between the beeps
     end
 end
 
 %% Save the whole workspace 
 % Performed every time a measurement is made   
 % ProcessedDir is not saved since a new one will be used each time the 
 % workspace is loaded
 % regexp: match regular expression (case sensitive)
 save(strcat(dataDir,pathSep,numberedFolderMeas,pathSep, ...
      'workspace.mat'),'-regexp','^(?!(ProcessedDir)$).');
  
elseif proc == 1 %  meas == 0 always in order to enter here
    
 %% Load previous measurement (whole workspace)
 % clearvars: these variables are used on the next load and then they are
 % replaced by the new ones that are loaded
 switch useLastMeas 
  case 0 % Doesn't load anything
     warning('useLastMeas=0 will not load things for proc=1')
  case 1 % Loads the last measurement
     clearvars -except dataDir pathSep numberedFolderMeas ProcessedDir; 
     load(strcat(dataDir,pathSep,numberedFolderMeas,pathSep,'workspace.mat'));
  case 2 % Loads a user-defined measurement
     clearvars -except dataDir pathSep measFoldName ProcessedDir; 
     load(strcat(dataDir,pathSep,measFoldName,pathSep,'workspace.mat'));
 end
 proc = 1; % In order to maintain its intended value (the program only enters to the 
           % switch if proc=1 and if by loading its value is changed, it gets restored to 1 here)
 meas = 0; % Same reason as above
 
end % End of measurements

%% Post-processing of the data and saving
if proc
    f_ProcessData(measfullpath,refmeasfullpath,ProcessedDir,pathSep, ...
    infoDelim,dataformat,cameraPlane,totalImgs,AiryFactor,metricSel, ...
    metricProfile,shiftCart);
  % Metric of the degree of extintion applied. Saves 1 plot per image
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% ACADEMIC PURPOSES: Zernike, simulation

%% Additional information
% abs(mask): should always be 1, meaning that it is normalized

%% Fourier transform of the mask
% Executed if desired on the parameters
% Performs the FFT of the mask and shows x and y profiles
f_ComputeMaskSpectrum(x,y,mask,abs_angFT,maskFTlog,FTmask);

%% Gradient of the mask
% Executed if desired on the parameters
% Shows the singularity clearly
AngMask = angle(mask);
f_ComputeMaskGradient(x,y,AngMask,gradMask);

%% Reconstruction of the mask with Zernike polynomials
% Executed if desired on the parameters
% Computes a wavefront reconstruction using Zernike's Polynomials and 
% calculates the function expansion coefficients
f_ZernikeReconstruction(z_ReconstrNumb,wrapMask,z_pupil,maskZernReconstr);

%% Simulation
% Executed if desired on the parameters
if simBool
 f_SimulateFreeSpace(x,y,Xpc,Ypc,rPC,mask,starAmplitude,planetAmplitude,...
 pixelSeparation,w1,w2,rPupilSize,showIin,showPupilin,showFPmag, ...
 logscale,showFPphas,showPhasout,showMagout,showIout)
end