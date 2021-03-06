%% Phase mask SLM displaying for the digital vortex coronagraphy technique
%% vortex images acquisition and processing.
% free-space simulation of two stars is also included optionally
%
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
%  -It is possible to load previous measurements of windows in a linux or
%  viceversa: with loadMeas = 3
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
%  MATLAB or vice versa, due to the incompatibility of the "pathSep"
%  -Errors or warnings may appear when you have for example a folder
%  "date-meas_2", but you don't have the folder "date-meas_1" and
%  "date-meas"; therefore, rename "date-meas_2" as "date-meas"
%
% Samuel Plazas(PA1/PA2/TG) - Juan Jose Cadavid(Master thesis) - 2018/2019



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PHASE MASK GENERATION WITH SLMs, VORTEX ACQUISITION VIA CAMERAS AND 
%%% PROCESSING WITH CORONAGRAPHIC METRICS
%% Parameters and initialization
%clc; clear; close all; % Initialization (SOMETIMES USED)
Parameters; % Adds to the algorithm all the needed parameters
% open Parameters; % Keep open always

% Default parameters for measuring: uncomment this in case you want to, but
% comment Parameters above
% ParametersMeasuring

%% Directories and add functions
% Adds initial functions
addpath(strcat('..',pathSep,toolsFldr,pathSep,filemanag)); 
[analysDir,toolsDir,dataDir,outDir] = ...
f_addDirectories(analysFldr,toolsFldr,dataFlrd,outFlrd,pathSep);
% Adds all the directories to use in the algorithm

%% Self-centering and spatial definition
if meas == 1 || (meas == 0 && proc == 0 ) || (measDebug == 1 && meas == 1)
    %% Self centering hardware initialization
    % This is here because the camera is needed for the self-centering
    % algorithm inside "f_DefineSpace"
    SLMcenterWisdom = strcat(dataDir,pathSep,'SLMwisdom.mat'); % ShiftCart 
    % is saved here for the self centering
    
    %% Self centering Lyot camera parameters
    % Always used for the self centering in this setup (of 2019)
    cameraLyot = 'DMK42BUC03';
    cameraIDLyot = 2;
    exposureLyot = 1/250; % Range: [1/1e4,1]
    fpsLyot = [];
    formatLyot = 'Y800 (1280x960)';
    cameraPlaneLyot =  'Lyot';
    PPLyot = 3.75; % Pixel pitch in [um/pixel]
    wisdomexist = (exist(SLMcenterWisdom,'file') == 2);
    if shiftMask == 2 && ~(wisdomexist)% For self-centering's execution
        %% Self centering camera preparation
        % Turns the camera on and create all the needed variables. Remember 
        % to leave the preview open
        [vidSelfCent,~,fpsLyot] = f_selectCamera(cameraIDLyot, ...
        cameraLyot,exposureLyot,fpsLyot,formatLyot);
    else % ~(measSimulated == 0 && meas == 1)
        vidSelfCent = []; % Empty, as a non-used input in some functions
    end
    
    %% Spatial definitions
    % ShiftCart will be used in "Measurement debugging"
    [rSize,x,y,Xslm,Yslm,rSLM,phiSLM,Xpc,Ypc,rPC,phiPC,monitorSize, ...
    shiftCart,mainLyotRadius] = f_DefineSpace(vidSelfCent,sSupport, ...
    sSize,shiftCart,shiftMask,PPLyot,pixSize,apRad,scrnIdx,circularMask,...
    z_pupil,coordType,MaxMask,SLMcenterWisdom,cameraLyot, ...
    cameraPlaneLyot,exposureLyot,formatLyot,fpsLyot,measSimulated,maskSel);
    % Defines the cartesian/polar coordinates, its sampling interval and
    % discretized angular part for gl
    %%% Coordinates selection:
    [X,Y,r,phi] = f_SelectCoordinates(Xslm,Yslm,rSLM,phiSLM,Xpc,Ypc,rPC,...
    phiPC,plotMask);
    % f_DefineSpace and f_SelectCoordinates were divided so that in
    % f_SelectCoordinates one just selects already-calculated variables
    
    %% Single phase mask selection and plot on the screen or on the SLM
    % Don't plot the mask if one will process or measure, but plot it for a
    % measurement debug (that requieres measSimulated = 0)
    if (proc == 1 || meas == 1) && ~(measSimulated == 0 && measDebug == 1)
        plotMask = 0;
    end
    [mask,wrapMask,maskFig,maskName] = f_PlotSelectedMask(X,Y,r,phi, ...
    phaseValues,tc,s,ph0,p,WsizeRatio,L,f_FR,bcst,period,T0,frkTyp, ...
    Aalpha,Angalp,Angbet,z_coeff,z_a,z_pupil,z_disp_wrap,z_plot,normMag,...
    binMask,binv,MaskPupil,rSize,monitorSize,scrnIdx,coordType,abs_ang, ...
    MaxMask,plotMask,maskSel);
    % Dependencies:
    % f_PlotSelectedMask -> f_SpiralMask (or any other) -> f_ProjectMask ->
    % f_MaskWrapCircDiscret -> f_discretizeMask -> f_wrapToRange
    %                       -> f_ScaleMatrixData
else
    [maskName] = f_getMaskName(maskSel);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% MEASUREMENT/PROCESSING BY AN AUTOMATED PARAMETER VARIATION
%% Folders and register creations on Data and Output    

[imgpartPath,measfullpath,ProcessedDir,ltcvect,lglvect,totalImgs, ...
numberedFolderMeas] = f_CreateFoldersRegisters(maskName,tcvect,glvect, ...
slm,cameraPlane,dataDir,outDir,pathSep,infoDelim,dirDelim,lastmeas, ...
meas,proc,loadMeas);

%% Measurement
if meas
  if measSimulated == 0 % When a real measurement will be
    % performed
    %% Camera preparation
    if measDebug
      imaqreset; % Disconnect and delete all image acquisition objects
      % Refresh image acquisition hardware by restoring the settings
    end
    % Turns the camera on and create all the needed variables. Remember to
    % leave the preview open
    [vidMeas,srcMeas,fps] = f_selectCamera(cameraID,camera,exposure, ...
                                           fps,format);
  else % ~(measSimulated == 0 && meas == 1)
    vidMeas = []; % Empty, just to be as a non-used input in some functions
    srcMeas = [];
  end
  
  %% Measurement debugging
  % Usefull for aligning the vortex and adjusting exposure parameters

  if measDebug == 1 && measSimulated == 0
    if plotMask ~= 2
      error('plotMask must be 2 in order to perform the debugging');
    end
    [~,circShiftX,circShiftY] = f_CaptureImage(vidMeas,dataDir, ...
    filename,imgformat,pathSep,infoDelim,dirDelim,snapsfldr, ...
    previewBool,loghist,camera,cameraPlane,exposure,format,fps, ...
    maskFig,plotMask);
    if circShiftX ~= 0 || circShiftY ~= 0 % doesn't change if it didn't
                                          % change at all
      [coorRelshiftX, coorRelshiftY] = f_circShiftToCart(circShiftX, ...
                                 circShiftY,monitorSize(1),monitorSize(2));
      [shiftX,shiftY] = f_calcHScoorToSgnCoor((coorRelshiftX + ... 
      monitorSize(1)/2)/monitorSize(1),(coorRelshiftY + ...
      monitorSize(2)/2)/monitorSize(2));
      shiftX = shiftX + shiftCart(1);
      shiftY = shiftY + shiftCart(2);
      shiftCart = [shiftX, shiftY];
      answer1 = input('ShiftCart was modified. Save changes? y/n: ','s');
      % returned always as a string  
      if strcmp(answer1,'y') % save
        save(SLMcenterWisdom,'shiftCart','-append')
      else
        disp('No changes were made.');
      end
    end
    % Takes a camera shot,shows a figure and saves it
    % Close the mask figure after debugging
    if ishandle(maskFig)
      close(maskFig);
    end
    
  %%% Get some hardware/software/tools info:
  % get(vid): displays the general parameters of the camera
  % getselectedsource(vid): an existing variable (src). Similar to get(vid)
  % imaqhwinfo(vid): displays the driver connection with MATLAB
  % imaqhwinfo: all the installed adaptors in the Image Acquisition Toolbox
  % disp(vid): displays acquisition information
  % imaqtool: toolbox for the camera
  % imaqreset: refresh image acquisition hardware by restoring the settings
  %            Disconnect and delete all image acquisition objects
  % propinfo(src): shows the allowed values that can be modified
    
    proc = 0; % Nothing should be processed when debugging
  else % Real measurement (measDebug == 0)
    
    %% Automated measurement
    % Performs the measurements and stores them:
    [refmeasfullpath] = f_AutomateMeasurement(srcMeas,vidMeas,Xslm,Yslm,...
    rSLM,phiSLM,Xpc,Ypc,rPC,phiPC,s,ph0,p,WsizeRatio,L,f_FR,bcst,period,...
    T0,frkTyp,Aalpha,Angalp,Angbet,z_coeff,z_a,z_pupil,z_disp_wrap, ...
    z_plot,normMag,binMask,binv,MaskPupil,rSize,monitorSize,scrnIdx, ...
    coordType,abs_ang,MaxMask,maskSel,ltcvect,lglvect,totalImgs, ...
    waitbeforemeas,imgpartPath,dataformat,imgformat,measfullpath, ...
    infoDelim,cameraPlane,tcvect,glvect,measSimulated,recordingDelay, ...
    whiteblackref,beepSound);
    
    %% Save the 'whole' workspace
    % Performed every time a measurement is made
    save(strcat(dataDir,pathSep,numberedFolderMeas,pathSep, ...
    'workspace.mat'),'measfullpath','refmeasfullpath','infoDelim', ...
    'dataformat','imgformat','cameraPlane','totalImgs','beepSound','L', ...
    'PP','mainLyotRadius','measSimulated','tcvect','glvect','slm','maskSel');
    % The saved variables include all the input ones of  f_ProcessData.m
    % minus the ones that are after "-except" in the "Load previous
    % measurement (whole workspace)" section
  end

%% Termination of the hardware 
f_releaseCamera(vidMeas); % Clear vid for the self-centering stage
  
elseif proc == 1 %  meas == 0 always in order to enter here (since here one
                 % will load a measurement)
  % Either it measures or loads a measurement
  
  %% Load previous measurement (whole workspace)
  % clearvars: these variables are used on the next load and then they are
  % replaced by the new ones that are loded
  switch loadMeas
   case 1 % Loads the last measurement
    clearvars -except dataDir pathSep numberedFolderMeas ProcessedDir ...
    metricSel metricProfile AiryFactor NA loadMeas;
    loadMeasPath = strcat(dataDir,pathSep,numberedFolderMeas,pathSep, ...
    'workspace.mat');
    if exist(loadMeasPath, 'file') == 2 % File exists
       load(loadMeasPath); 
    else % File does not exist
      error('Could not load the specified previous measurement')  
    end
    
   case 2 % Loads a manually-put folder name
    clearvars -except dataDir pathSep measFoldName ProcessedDir ...
    metricSel metricProfile AiryFactor NA loadMeas;
    loadMeasPath = strcat(dataDir,pathSep,measFoldName,pathSep, ...
                          'workspace.mat');
    if exist(loadMeasPath, 'file') == 2 % File exists
       load(loadMeasPath); 
    else % File does not exist.
      error('Could not load the specified previous measurement') ;
    end
    
   case 3 % Loads a user-defined measurement
    % Selects only the folder from the full path
    clearvars -except dataDir pathSep ProcessedDir ...
    metricSel metricProfile AiryFactor NA loadMeas ;
    %%%  Select a custom measurement folder
    dataFlrd = 'Data';
    usermeasFoldName = uigetdir; % The user will select the directory
    [folderRoot,usermeasFoldName] = fileparts(usermeasFoldName);
    newDirBase = strcat(folderRoot,pathSep,usermeasFoldName,pathSep);
    loadMeasPath =  strcat(newDirBase,'workspace.mat');
    
    if exist(loadMeasPath, 'file') == 2 % File exists
       load(loadMeasPath); 
    else % File does not exist.
      error('Could not load the specified previous measurement')  
    end
    measfullpath =  strcat(newDirBase,'PSF-allmeas.mat'); % Update new dir
    refmeasfullpath = strrep(strcat(folderRoot,pathSep,refmeasfullpath(...
    strfind(refmeasfullpath,dataFlrd)+length(dataFlrd)+1:end)),'\', ...
    pathSep);
    
   otherwise % error
    errror('Invalid value for "loadMeas"')
  end
  
  %%% Keep some input variables defaul:
  proc = 1; % In order to maintain its intended value (the program only
  % enters to the switch if proc=1 and if by loading its value is
  % changed, it gets restored to 1 here). Ref: proc = 1
  meas = 0; % Same reason as above. Ref: meas = 0
  
end % End of measurements

%% Post-processing of the data and saving
% It is possible to measure and then process
% It is also possible to load a measurement and process it
if proc == 1% && measSim == 0
  if ~exist('refmeasfullpath','var')
     error(['"refmeasfullpath" does not exist, no previous' ...
            ' measurement was found']);
  end
  % Metric of the degree of extintion applied
  [x,y,radialIntensityRef,radialIntensityMeas,cartcoord,refEEFcurve, ...
  refNormIntensity,measEEFcurves,measNormIntensity,logSNR, ...
  attenuationRatio,EEFattenuationRatio,refIntensityGradient, ...
  measIntensityGradient] = f_ProcessData(measfullpath,refmeasfullpath, ...
  ProcessedDir,dataDir,pathSep,totalImgs,metricProfile, beepSound,L,NA, ...
  PP,measSimulated,glvect,tcvect);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% ACADEMIC PURPOSES: Zernike, simulation
if proc == 0
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
  f_ZernikeReconstruction(z_ReconstrNumb,wrapMask,z_pupil, ...
                          maskZernReconstr);

  %% Simulation
  % Executed if desired on the parameters
  if simBool
    f_SimulateFreeSpace(x,y,Xpc,Ypc,rPC,mask,starAmplitude, ...
    planetAmplitude,pixelSeparation,w1,w2,rPupilSize,showIin, ...
    showPupilin,showFPmag,logscale,showFPphas,showPhasout,showMagout, ...
    showIout)
  end
end