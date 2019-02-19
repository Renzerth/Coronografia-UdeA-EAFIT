%%  Phase mask SLMs display for the digital vortex coronagraphy and 
%%% free-space simulation of two stars
% Inputs:
% See "Paramters.m"
%
% Outputs:
%  Plots on the PC or on the SLMs
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
% Folders:
%  Analysis: principal scripts
%  Data: the inputs of the algorithm are the acquired vortex images
%  Data -> DatalogDir: specific measurement folder 
%  Output: processed images or plots
%  Output -> ProcessedDir: specific processed images folder
%  Tools: functions used in the program
%
% Bugs:
%  When the ticks of -pi, ..., pi are shown, sometimes the may move
%  slightly from the real value
%  The algorithm shall only be executed inside the analysis folder
%  Not tested yet in Linux
%  Simulations are not fully accurate with the calculations
%  Zernike polynomials are not fully wrapped and the "a" constant has not
%  been fully understood
%
% Samuel Plazas(PA1/PA2/TG) - Juan Jose Cadavid(Master thesis) - 2018/2019

%%%%%%%%%%%%%%%%%%%%%%%%%%% PHASE MASK GENERATION ON THE SLM's 
%% Parameters and initialization
clc; clear; close all; % Initialization
Parameters; % Adds to the algorithm all the needed parameters
% open Parameters; % Keep open always

%% Directories and add functions
addpath('sub_scripts'); % Adds all the sub programs of the algorithm
addDirectories; % Adds all the directories to use in the algorithm

%% Spatial definitions
spatialDefinitions; % Defines the cartesian/polar coordinates, its 
                    % sampling interval and discretized angular part for gl
                
%% Plot one mask for tests: phase mask selection
PhaseMaskSel; % Selects the type of phase mask.
              % The variables "mask" and "maskName" are outputs here

%% Phase mask plot on the screen or on the SLM
%%% Plot phase mask on the Fourier plane of the vortex coronagraph
f_fig_maskSLM(x,y,r,mask,gl,glphi,mingl,maxgl,levShft,abs_ang,binMask,monitorSize,plotMask);




%%%%%%%%%%%%%%%%%%%%%%% MEASUREMENTS BY AN AUTOMATED PARAMETER VARIATION
if meas == 1
%% Folders and register creations on Data and Output    
FoldersRegistersCreation;

%% Hardware initialization
% HardwareInit; % Future script 
                % Turns the camera on and create all the needed variables
                % Remember to leave the preview open
[vid,src] = f_selectCamera(camera,exposure,format);
% USE: ???
% - vid.FramesPerTrigger = 1; % Default frames to capture per trigger
% - 

%% Measurement debugging
% Usefull for aligning the vortex and adjusting exposure parameters
if measDebug == 1
    f_ImageCapture(vid,dataDir,filename,imgformat); % Takes a camera shot,
                                                    % shows a figure and 
                                                    % saves it
                                          
    Frame = f_GetFrame(vid); 
else
    %% Reference measurement
    % Still not sure if needed: null tc beam or a high tc beam(long radius)

    %% Automated measurement
    AutomatMeasure; % Future script

    %% Post-processing of the data
    DataProcessing; % Metric of the degree of extintion applied

    %% Save data
    % SaveData; % Future script
end

%% Termination
% Terminate_settings; % Future script % Clears variables and closes all; deactivates camera
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
% abs(mask) should always be 1, meaning that it is normalized; try by yourself

%% Fourier transform of the mask
if FTmask == 1
    maskFT; % Performs the FFT of the mask and shows x and y profiles
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