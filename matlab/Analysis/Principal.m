%% Phase Mask SLMs Display and Free-Space Simulation of two stars
% Inputs:
% See "Paramters.m"
%
% Outputs:
%  Plots on the PC or on the SLMs
%  Processed vortex images
%
% Notes:
%  Units: a.u (arbitrary units) and cm for lengths, radians for angles and
%  um for wavelengths
%  Spiral phase masks are inside Laguerre-Gauss
%  The variable mask is complex and is wrappped: mask = exp(i*mask)
%  All wrappped phases are shown on [-pi,pi] 
%  Always execute the program whenver you are exactly inside its folder
%
% Samuel Plazas Escudero - Juan Jose Cadavid - 2018/2019 - PA1/PA2/TG

%%%%%%%%%%%%%%%%%%%% PHASE MASK GENERATION ON THE SLM's %%%%%%%%%%%%%%%%%%%

%% Parameters and initialization
clc; clear; close all; % Initialization
Parameters; % Adds to the algorithm all the needed parameters
% open Parameters; % Keep open always

%% Directories and add functions
analysDir = pwd; cd ..; % Store script directory
cd(toolsFldr); toolsDir = pwd; cd ..; % Store function directory
cd(dataFlrd); dataDir = pwd;  cd ..; % Store data directory
addpath(genpath(toolsDir)); cd(analysDir); % Add all folders in functions
% restore back default paths, type: restoredefaultpath

%% Spatial definitions
sSize = 2^k-1; % Number of samples; odd number so that vortex gets
               % centered (spatial size); Spatial size. ref: 2^k-1
SpatialSupport = SpatialSupport/2; % Half support of the SLM window in cm
spaceVector = -SpatialSupport:2*SpatialSupport/(sSize-1):SpatialSupport; % Symmetric space
[X,Y] = meshgrid(spaceVector); % A symmetric grid: 2D Cartesian coordinates
[phi,r] = cart2pol(X,Y); % Polar coordinates
x = spaceVector; % Cartesian x-vector
y = x; % Cartesian y-vector

%% Plot one mask for tests: phase mask selection
PhaseMaskSel; % Selects the type of phase mask.
              % The variables "mask" and "maskName" are outputs here

%% Phase mask plot on the screen or on the SLM
%%% Plot phase mask on the Fourier plane of the vortex coronagraph
f_fig_maskSLM(x,y,r,mask,m,n,a,b,gl,abs_ang,binMask,plotMask);

%%%%%%%%%%%%%%%%% Automated Measurements varying parameters%%%%%%%%%%%%%%%%
if meas == 1
%% Folder creation
ltcvect = length(tcvect); % Length of the tc vector
lglvect = length(glvect); % Length of the gl vector
strDate = date; % Today's date is retrieved from the local machine
MeasSize = [maskName '_mask_tcs_' num2str(ltcvect) '_gls_' num2str(lglvect)]; % Datalog with the number of measurements for tc's and gl's
Datalogfldr = [date '_' MeasSize]; % Folder name
cd(dataDir);
ax = exist(Datalogfldr, 'dir'); % 7 if folder exists, 0 if not
if ax ~= 7 % Create a folder if it doesn't exist
    mkdir(Datalogfldr);
end
cd(analysDir);

%% Hardware initialization
HardwareInit; % Turns the camera on and create all the needed vars
              % Remember to leave the preview open

%% Measurement debugging
% f_CameraShot(); % Takes a photo, shows a figure and saves it as shot.png
% Usefull for aligning the vortex and adjusting exposure parameters

%% Reference measurement
% tc = 0 beam or a high tc beam (long radius)

%% Automated measurement
% AutomatMeasure
%addpath(PhaseMaskSel);
cd(dataDir); % Go to data Directory: input for the program is the acquired 
             % images of the vortices
showM = 0; % Don't show a fig in "PhaseMaskSel"
wait = 2; % Waits 2 second between each mask to be shown

for i = 1:length(tcvect)
  for j = 1: length(glvect)
    tc = tcvect(i); % Specific tc for this iteration
    gl = glvect(j); % Specific gl for this iteration
    cd(analysDir); % Moves to the Analysis directory
    PhaseMaskSel; % Selects a phase mask to display
    f_fig_maskSLM(x,y,r,mask,m,n,a,b,gl,abs_ang,binMask,plotMask);
    pause(wait); % Displays the mask for "wait" seconds   
    MeasInfoEach = ['tc_' num2str(tcvect(i)) '_gl_' num2str(glvect(j))];
    cd(dataDir); cd(Datalogfldr); % Goes to the data directory and
                                  % enters to the specific measurement
                                  % folder
    saveas(gcf,[MeasInfoEach '.png']); % Saves the last shown figure
                                    % plotMask should be different from 0
    
    % Right now, the masks are being saved, but later the images should
    % be saved as an input of the algorithm to be processed
    % f_CameraShot();
     a=1;   
  end
end

cd(analysDir); % Go to script Directory

% Maybes:
% Calculate the energy of the images
% Generalize the way the fors are executed
% Randomize the order of the tcvect and glvect, commenting the premises of
% the experimental design on PA2

%% Post-processing of the data
% DataProcessing


%% Save data
% SaveData

%% Termination
% Terminate_settings; % Clears variables and closes all; deactivates camera
end
















%%%%%%%%%%%%%%%%%%%% NOT USED BUT FOR ACADEMIC PURPOSES %%%%%%%%%%%%%%%%%%%
% abs(mask) should always be 1, meaning that is normalized; try by yourself

%% Optional FT
if FTmask == 1
    FFT2D = @(s) ifftshift((fft2(fftshift(s)))); % 2D Fourier Transform
    mask = FFT2D(mask);
    f_fig_maskSLM(x,y,r,mask,m,n,a,b,gl,abs_ang,binMask,plotMask);
end

%% Reconstruction of the mask with Zernike polynomials
if maskZernReconstr == 1
    f_Zernike_Reconstruction(14,angle(mask),1);
end

%% Simulation
if sim == 1
Simulation; % Executed if desired on the parameters
end