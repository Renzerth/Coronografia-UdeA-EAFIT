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
%  Units: a.u (arbitrary units) and cm for lengths, radians for angles and
%  um for wavelengths
%  Helicoidal phase masks are inside the Laguerre-Gauss beams category
%  The variable mask is complex and is wrappped: mask = exp(i*mask)
%  All wrappped phases are shown on [-pi,pi] 
%  Always execute the program whenver you are exactly inside its folder
%
% Bugs:
%  -When the ticks of -pi, ..., pi are shown, sometimes the may move
%  slightly from the real value
%  -The algorithm shall only be executed inside the analysis folder
%  - Algorithm not tested yet in Linux but in Windows 7 and 10
%  -Simulations are not fully accurate with the calculations
%  -Zernike polynomials are not fully wrapped and the "a" constant has not
%  been fully understood
%  -If a camera is not recognized or even the screen index are not acting
%  coherently, restart MATLAB
%  -For high resolution monitors, sometimes the figure bar may appear but
%  this presents no problem for the mask projection
%  
%
% Samuel Plazas(PA1/PA2/TG) - Juan Jose Cadavid(Master thesis) - 2018/2019

%%%%%%%%%%%%%%%%%%%%%%%%%%% PHASE MASK GENERATION ON THE SLM's 
%% Parameters and initialization
clc; clear; close all; % Initialization
Parameters; % Adds to the algorithm all the needed parameters
% open Parameters; % Keep open always

%% Directories and add functions
% Adds initial functions
addpath(strcat('..',pathSep,toolsFldr,pathSep,filemanag)); 
[analysDir,toolsDir,dataDir,outDir] = ...
f_addDirectories(analysFldr,toolsFldr,dataFlrd,outFlrd);
% Adds all the directories to use in the algorithm

%% Spatial definitions
[Xslm,Yslm,rSLM,phiSLM,Xpc,Ypc,rPC,phiPC,sSize,monitorSize] = ...
f_DefineSpace(spaceSupport,shiftCart,pixSize,scrnIdx,circularMask, ...
shiftBool,coordType);
% Defines the cartesian/polar coordinates, its sampling interval and 
% discretized angular part for gl
                
%% Phase mask selection and plot on the screen or on the SLM
[mask,maskName] = f_PlotSelectedMask(Xslm,Yslm,rSLM,phiSLM,Xpc,Ypc,rPC,phiPC,gl, ...
glphi,mingl,maxgl,levShft,tc,s,ph0,p,W,binv,norm, L,f_FR,bcst,z_coeff, ...
a,frac,pupil,sSize,disp_wrap,plot_z,binMask,monitorSize,scrnIdx, ...
abs_ang,plotMask,maskSel);   
% Dependencies:
% f_PlotSelectedMask -> f_SpiralMask (or any other) -> f_ProjectMask -> 
% f_MaskWrapCircDiscret -> (f_discretizeMask & f_ScaleMatrixData)          

                           
                           
%%%%%%%%%%%%%%%%%%%%%%% MEASUREMENTS BY AN AUTOMATED PARAMETER VARIATION
if meas == 1
 %% Folders and register creations on Data and Output    
 CreateFoldersRegisters;

 %% Hardware initialization
 if measSimulated == 0 % When a real measurement will be performed
  % InitializeHardware;
               % Turns the camera on and create all the needed 
                     % variables. Remember to leave the preview open
  [vid,src] = f_selectCamera(camera,exposure,format);
 end

 %% Measurement debugging
 % Usefull for aligning the vortex and adjusting exposure parameters
 if measDebug == 1
  SingleFrame = f_CaptureImage(vid,dataDir,filename,imgformat,pathSep, ...
                               snapsfldr); 
  % Takes a camera shot,shows a figure and saves it   
  figure; imhist(SingleFrame); % Shows a histogram of the snapshot

  % Get some hardware/software/tools info:
  % get(vid): displays the general parameters of the camera
  % getselectedsource(vid): an existing variable (src). Similar to get(vid)
  % imaqhwinfo(vid): displays the driver connection with MATLAB
  % imaqhwinfo: all the installed adaptors in the Image Acquisition Toolbox
  % disp(vid): displays acquisition information
  % imaqtool: toolbox for the camera
  % imaqreset: refresh image acquisition hardware by restoring the settings
 else
  %% Reference measurement
  % Still not sure if needed: null tc beam or a high tc beam(long radius)

  %% Automated measurement
  AutomateMeasurement; % Performs measurements and stores them

  %% Post-processing of the data and saving
  %ProcessData; % Metric of the degree of extintion applied
                % Saves plot(s) of the applied metrics
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
end % End of measurements




%%%%%%%%%%%%%%%%%%%%%%% ACADEMIC PURPOSES: Zernike, simulation
% abs(mask): should always be 1, meaning that it is normalized

%% Fourier transform of the mask
if FTmask == 1
    maskSpectrum; % Performs the FFT of the mask and shows x and y profiles
end

%% Gradient of the mask
% Shows the singularity clearly
if gradMask == 1
    [xg,yg] = gradient(angle(mask));
    figure; contour(x,y,angle(mask)); hold on;
    quiver(x,y,xg,yg); title('Gradient of the mask'); hold off
    figure; imagesc(x,y,xg); 
    colormap(hot); title('X-profile gradient of the mask');
    figure; imagesc(x,y,yg); 
    colormap(hot); title('Y-profile gradient of the mask');
end

%% Reconstruction of the mask with Zernike polynomials
if maskZernReconstr == 1
    f_Zernike_Reconstruction(14,angle(mask),1);
end

%% Simulation
if sim == 1
Simulation; % Executed if desired on the parameters
end