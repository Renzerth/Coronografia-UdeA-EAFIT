%% Phase Mask SLMs Display and Free-Space Simulation of two stars
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
%  Data -> Datalogdir: specific measurement folder 
%  Output: processed images or plots
%  Tools: functions used in the program
%
% Bugs:
%  So far no
%
% Samuel Plazas (PA1/PA2/TG) - Juan Jose Cadavid(Master thesis) - 2018/2019

%%%%%%%%%%%%%%%%%%%%%%%%%%% PHASE MASK GENERATION ON THE SLM's 
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




%%%%%%%%%%%%%%%%%%%%%%% MEASUREMENTS BY AN AUTOMATED PARAMETER VARIATION
if meas == 1
%% Measurement folder creation (Datalog)
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
Datalogdir = [dataDir '/' Datalogfldr]; % Specific measurement folder
cd(analysDir);

%% Hardware initialization
% HardwareInit; % Future script % Turns the camera on and create all the needed vars
                % Remember to leave the preview open
%[vid,src] = f_selectCamera(camera,exposure,format);
% Use vid.FramesPerTrigger = 1; ??


%% Measurement debugging
% f_ImageCapture(vid,dataDir,filename);
% Frame = f_GetFrame(vid);

% f_CameraShot(); % Future script % Takes a photo, shows a figure and saves it as shot.png
% Usefull for aligning the vortex and adjusting exposure parameters

%% Reference measurement
% Null tc beam or a high tc beam (long radius)

%% Automated measurement
% AutomatMeasure; % Future script
showM = 0; % Don't show a fig in "PhaseMaskSel"
recordingDelay = 2; % Waits 2 second between each mask to be shown
totalImgs = ltcvect*lglvect; % Number of images to be taken
expImgs = cell(1,totalImgs); % Cell with the experimental images
MeasInfo = expImgs; % Same initialization as expImgs
idxgral = 1; % General index that will be on [1,totalImgs]
% fileFormat = '.bmp'; % OLD: save each image

for idxtc = 1:ltcvect 
  for idxgl = 1:lglvect
    tc = tcvect(idxtc); % Specific tc for this iteration
    gl = glvect(idxgl); % Specific gl for this iteration
    cd(analysDir); % Moves to the Analysis directory
    PhaseMaskSel; % Selects a phase mask to display
    f_fig_maskSLM(x,y,r,mask,m,n,a,b,gl,abs_ang,binMask,plotMask);
    pause(recordingDelay); % Displays the mask for "recordingDelay" seconds   
    tcstr = ['tc_' num2str(tcvect(idxtc))]; 
    glstr = ['_gl_' num2str(glvect(idxgl))];
    MeasInfo{idxgral} = [tcstr glstr]; % Dataname for each experimental data
    cd(Datalogdir); % Goes to the data log directory (specific measurement
                             % folder)
    % snap = getsnapshot(vid); % Real measurements
    wrappedMask = f_circularPupil_maskAngle(r,mask,binMask); 
    snap = wrappedMask; % "Simulated" measurements (the mask is saved)
    expImgs{idxgral} = snap;
    
    % save(filename,variables,'-append')
    A=expImgs{idxgral};
    save(MeasInfo{idxgral},'A');
    % imwrite(expImgs{idxgral},[MeasInfo{idxgral} fileFormat]); % OLD: save each image % Saves the last shown figure
    % Here, the masks are saved, but later the images should
    % be saved as an input of the algorithm to be processed
    % f_CameraShot();
    idxgral = idxgral + 1; % The general index increases   
  end
end

%save(MeasInfo{idxgral},'expImgs')

cd(analysDir); % Go to script Directory

%% Post-processing of the data
% DataProcessing; % Future script


%% Save data
% SaveData; % Future script

%% Termination
% Terminate_settings; % Future script % Clears variables and closes all; deactivates camera
% delete(vid); % Clean up the camera

%% End notification
if beepSound == 1
    for beepTimes = 1:3
        beep();
        pause(0.2);
    end
end

end




%%%%%%%%%%%%%%%%%%%%%%% ACADEMIC PURPOSES: Zernike, simulation
% abs(mask) should always be 1, meaning that it is normalized; try by yourself

%% Optional FT
if FTmask == 1
    FFT2D = @(s) ifftshift((fft2(fftshift(s)))); % 2D Fourier Transform
    mask = FFT2D(mask); % FT of the mask (not wrapped)
    mask = abs(mask); % Magnitude of the FT
    mask = 20*log(mask.^2); % Magnitude squared of the FT in log scale
    abs_ang_FT = 1; % Magnitude is always plotted
    f_fig_maskSLM(x,y,r,mask,m,n,a,b,gl,abs_ang_FT,binMask,plotMask);
    
    [maxX, maxY] = size(mask);
    midX = ceil((maxX+1)/2);
    midY = ceil((maxY+1)/2);  
    
    figure;
    f = improfile(mask,[1,1023],[512,512]);
    plot(x,f); 
    title('Horizontal profile of abs squared FFT');
    
    figure;
    g = improfile(mask,[512,512],[1,1023]);
    plot(x,g); 
    title('Vertical profile of abs squared FFT');
end

%% Reconstruction of the mask with Zernike polynomials
if maskZernReconstr == 1
    f_Zernike_Reconstruction(14,angle(mask),1);
end

%% Simulation
if sim == 1
Simulation; % Executed if desired on the parameters
end