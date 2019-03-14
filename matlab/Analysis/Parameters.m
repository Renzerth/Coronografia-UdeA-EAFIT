%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% PART 1: GENERAL ADJUSTMENTS %%%%%%%%%%%%%%%%%%%%%%%
%% Algorithm sections
meas = 1; % Measure: yes (1) or no (0)
%%% For meas = 1:
    % Note: measDebug will be 0 if measSimulated = 1
    % If meas = 1 -> all the figures will be closed before starting it
    measDebug = 0; % Debugging before actually measuring. Displays the 
                   % default phase mask and shots a photo with the camera
                   % Works if  measSimulated = 0
    measSimulated = 1; % Saves the mask and does not involve the cameras: 
                       % yes (1) or no (0)
    beepSound = 1; % Beep sound when the measurement finishes.
proc = 0; % Processes the data

%% General algorithm parameters: coordinates, plots, screens and mask type
precision = 3; % Precision of displayed results: significative digits (3)
abs_ang = 2; % Custom(0)[str has to be defined for this case], magnitude
             % (1) or phase (2) plot. It doesn't apply for Zernike and LG +
             % Zernike: instead use z_disp_wrap for phase wrapping or not.
             % This is done since the aberration effects observed in the 
             % phase, are noticed in the amplitude if the field is 
             % propagated. Consider using simBool = 1
             % abs_ang = 0 is not valid for FTmask = 1
maskSel = 1; % Phase mask selection:
             % 0: Helicoidal mask: SPP or DSPP depending on gl
             % 1: Laguerre-Gauss beams: amplitude or phase
             % 2: VPL: Vortex Producing Lens = Helicoidal + Fresnel lens
             % 3: Elliptic Gaussian beam phase mask
             % 4: Fork phase mask
             % ---- NOT USED:
             % 5: Zernike (aberrations)
             % 6: Laguerre-Gauss + Zernike
             % 7: Hermite-Gauss beams NOT DONE
             % 8: Mutliple vortices NOT DONE
             % 9: Sum of spiral phase masks NOT DONE
             % 10: Gerchberg-Saxton NOT DONE
             % otherwise: Unitary
plotMask = 1; % Allows to plot the final mask, as it can be a combination 
              % of the previous ones
              % 0: no plot;
              % 1: on the screen
              % 2: on the SLM
              % 3: on the screen but surface-plot type
        
                     
                     
              
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%              
%%%%%%%%%%%%%%%%%%%%%%% PART 2: HARDWARE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% scrnIdx: screen number selector. In [1,N] with N the # of screens
% Windows 7 PC used in 2019 (according to):
%  -Principal screen: MATLAB scrnIdx(1); Windows(2); AnyDesk(2)
%  -Pluto screen: MATLAB scrnIdx(3); Windows(1); Anydesk(1)
%  -LC2002 screen: MATLAB scrnIdx(2); Windows(3); Anydesk(0)
slm = 'LC2002'; % 'Pluto' (reflection); 'LC2002' (transmission); 'No-SLM'
switch slm
  case 'Pluto'
    %% SLM parameters (reflection)
    % spaceSupport = 1; % Unitary space: spaceVector = -1:2/(Ssize-1):1;
    sSupport = min([1.536 0.864]); % Size of the SLM window in cm:
                                   % 1.536cm x 0.864cm
    maxNumPix = max([1920 1080]); % Maximum number of pixels on the SLM 
                                  % (either horizontal or vertical); SLM's 
                                  % resolution in pixels: 1920 x 1080 
    pixSize = 8; % SLM pixel's size in um
    scrnIdx = 3; % Screen number selector. Default: 3
    
  case 'LC2002'
    %% SLM parameters (transmision)
    sSupport = min([2.66 2.00]); % Same as the reflection SLM
    maxNumPix = max([800 600]); % Same as the reflection SLM
    pixSize = 32; % Same as the reflection SLM in um
    scrnIdx = 2; % Screen number selector. Default: 2
    
  case 'No-SLM'
   sSupport = 1; % A unitary space can be created when coordType=1
   maxNumPix = max([1920 1080]); % As the Pluto SLM
   pixSize = 8; % As the Pluto SLM
   scrnIdx = 1; % PC screem
  otherwise
    warning('Please select an SLM');
end 

%% SLM positionining calibration, coordinates and type of truncation
MaskPupil = 1; % Applies a pupil truncation to the mask: (0): no; (1): yes
% Won't work for maskSel = 5 or 6 (Zernike), as it has z_pupil
coordType = 1; % Type of calculation of the spatial coordinates. def: 2 
% 1: size defined by the user, space support defined by the SLM to use
% 2: size defined by the resolution of the selected screen    
%%%% For coordType = 1 (user custom-sized):
    k = 9; % Bits for grey levels; 2^k is the resolution (size of x and y)
           % Default: 10. Size is calculated as 2^k - 1 or 2^k in sSize
           % Usually try maximum k = 11
           % Only works when coordType = 1
    sSize = 2^k-1;  % Spatial size: number of samples; odd number so that 
                    % the vortex gets centered. ref: 2^k-1      
    % sSize = 512; % Instead of 2^k-1
    MaxMask = 1; % Defines if the mask should be maximized
    % 0: custom-size mask that depends on the variable sSize   
    % 1: maximizes the mask for coordType = 1
    % 2: maximized mask but keeping its rectangular fashion. MaxMask = 2 is
    %    analog to having circularMask = 1 for coordType = 1
%%% For plotMask=2 (SLM plotting):
    circularMask = 1; % Works either on coordType = 2 or when MaxMask = 1
     % 0: The mask presents an elliptical form when in the full screen
     % 1: The mask presents a circular form when in the full screen
     % On both cases full screen means that plotMask = 2
    shiftMask = 1; % Shift for all masks
     % 0: shift deactivated [for exporting masks]
     % 1: shift activated [SLM displaying]
     % 2: self-centering algorithm
    shiftCart = [31.5,-1.8]; % shiftCart: [yshift,xshift], works when 
                       % shiftMask = 1. Percentages of movement of the
                       % total size of the mask (cartesian coordinates 
                       % convention). Ranges per shift: [0,100] (percentge)  
                          
   % LC2002: [-3,0.1]
   % Pluto: [31.5,-1.8]
              
%% Camera selection and parameters
camera = 'DMK23U445';
% Exposure: analog parameter
% Format: 'Y800 (1280x960)' [best]; 'RGB24 (1024x768)' [another option]

switch camera
  case 'DMK23U445' % PSF plane % CMOS
    exposure = 1/1e3; % Range: [,]
    format = 'Y800 (1280x960)'; 
    cameraPlane =  'PSF';
    PP = 3.75; % Pixel pitch in um
    
  case 'DMK42BUC03' % Lyot plane % CCD
    exposure = 1/1e3; % Range: [1/1e4,1]
    format = 'Y800 (1280x960)';
    cameraPlane =  'Lyot';
    PP = 3.75; % Pixel pitch in um
    
  case 'DMK41BU02.H' % not used here % CCD
    exposure = 1/1e3; % Range: [,]
    format = 'Y800 (1280x960)';    
    cameraPlane = 'notusedhere';
    PP = 4.65; % Pixel pitch in um
end

% For 'DMK42BUC03' [Delete]:
% -Default: 0.0556 or 1.282
% -PH50um: 0.0030 -- Y0.0256 || PH25um: 0.0083 -- Y0.0556 ||
% -PH15um: --Y0.1429 || 0.0227 || PH10um: 0.0435 -- Y0.200 || 
% -PH5um: 0.363 -- Y3.099

%% Image capture
filename = 'test'; % Name of the capture one wants to take
imgformat = '.png'; % Format with period. mat, bmp, png, jpg
                    % This format doesn't apply for the measurements

                    
                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          
%%%%%%%%%%%%%%%%%%%%% PART 3: PHASE MASKS PARAMETERS %%%%%%%%%%%%%%%%%%%%%%
%% Parameters: Laguerre-Gauss, spiral phase mask and general masks
L = 0.6328; % Laser wavelength [um]. Used in Zernike and VPL masks
tc = 1; % Topological charge (integer bigger or equal to one)
        % tc = Azimuthal index m for LG. Fractional tc result on phase
        % patterns of Hermite-Gauss (maybe just a coincidence)
s = 1; % Sign of mask (+1 or -1); reverses the imprinted OAM 
ph0 = pi; % Initial phase of the angle [radians]; reference +pi from
         % normal zero of trig circle and same rotation convention.
         % This corresponds to a normal rotation of the mask for stethic
         % reasons and shouldn't affect the results. Only affects if the
         % vortex is no fully centered
% For abs_ang = 2:
binMask = 0; % Binarizes the mask w.r.t the max/min of the phase (boolean)
binv = 0; % Binary inversion of the mask: yes(1); no(0). Only applies when 
          % binMask=1. It is usefull to be applied for odd p's on LG beams
% For abs_ang = 1:
normMag = 0; % Normalize magnitude. yes(1); no(0). 
          
%% Parameters: Laguerre-Gauss
p = 1; % Number of radial nodes. If p=0, normal helicoid masks are obtained
       % If they are used and tc=0(m=0); binary masks are obtained
       % Even p; rings are ones. Odd p; rings are zeroes. Use mask = mask'
WsizeRatio = 100; % Width of the modes; for LG; ref: [0,100] (percentage  
                % that is more understandable when p=1)
                % When smaller than 3, the size of the mask decreases
         
%% Parameters: VPL Phase mask, 
% f_FR: Fresnel lens focal distance or diffractive lens phase focal length
f_FR = 1.3; % In m [the bigger, the more plane the phase is]
%%% Fixed parameters for the minimum focal length criteria:
AreaSLM = maxNumPix*pixSize^2; % SLM's area of the longest dimension [um]
minf_FRum = AreaSLM/L; % Criterium to determine the MINIMUM f_FR [um]
scaleFactor = 1e-6; % um to m. Constant factor
minf_FR = minf_FRum*scaleFactor; % um to m
if f_FR < minf_FR % f cannot be smaller than the criterium of the smallest
  error(['VPL criterium not fulfilled, establish f_FR bigger than ' num2str(minf_FR)]);
end
% From: 2_edgar_2015_Generation_Optical_Vortices_Binary_Vortex_Lenses.pdf
% minf_FR is used in meters in this paper and uses reasonable scales:
% 1_edgar_2013_High-quality optical vortex-beam generation_E-Rueda_OL.pdf
% The scale for working with the SLMs is meters

%% Parameters: Elliptic Gaussian Vortex
bcst = 0.3; % Ellipticity. cy/cx = 1/alpha. Ref: .1, .2, .4, .6, .8 and 1
            % Out of theory: when beta >> 1, the mask tends to be binary
            % When beta ~ 0, the mask tends to be trinary

%% Parameters: Fork Phase
frkTyp = 1; % 1: smooth transition (amplitude); 
            % 2: phase jump transition (phase)
            % 1 = 2 when both are binarized
period = 0.1; % Period of the grating (fringe spacing). Ref: 0.1
              % Twice its inverse is the number of line dislocations
              % So the frequency is f = 2/period
              % Recommended: L or L/2 (close to a real diffraction grating)
              % L: Laser wavelength [um]
%%% Smooth transition (amplitude filter)
T0 = 1; % Const. absorption coeff of the hologram; only affects amplitude. 
Aalpha = pi; % Amplitude of the phase modulation. Ref: pi
% Its value was found experimentally and allows the mask to have the
% desired range: [-pi,pi] (otherwise it starts to break for bigger values 
% or doesn't sample all the phase range for lower values)
% *Symmetric fork: if change ph0 is an odd multiple of pi/2. The line
% dislocation should align with the diffraction grating
% *Fork facing upwards: adjust either the sign of period or of s to be (-)
% *Binary inversion: changes for every odd multiple of pi/2 applied to ph0
%  or by changing the sign of alpha
% *Ramifications of the fork: equals to tc
% The charge of the vortex can be determined by counting the number of 
% forks, or subtracting one from the number of prongs.

%%% Phase jump transition (phase filter)
Angalp = pi/2; % Diffraction angle of horizontal direction (x) [radians]
Angbet = 0; % Diffraction angle of vertical direction (y) [radians]
% Here, more than pi/2 seems not to work very well
% Range: [-pi/2,pi/2]
% The "*" items remarks in the 'Smooth transition' comments also apply here

%% Parameters: Zernike polynomials with Noll's convention
%%%% For maskSel = 5 or 6:
z_coeff = [4]; % Zernike coeffient vector (see f_ZernikeMask.m)
z_a = 2.5; % Arbitrary constant; the bigger, the more intense; ref: a=2.5
z_pupil = 1; % Pupil relative size: [0,1]; like a percentage
z_disp_wrap = 1; % (0): Original; (1): wrapped mask on [-pi,pi] 
z_plot = 0; % plot with Zernike builder: yes(1); no(0)
% L and gl are also used with Zernike
% Test orthogonality:
% set: z_coeff = [0 0 0 0 0.5 0 0 0 0 0 0 ] and open f_ZernikePolynomials
% and uncomment the "Kronecker Delta Plotting " section
% For coordType = 1, a low value of k (for example less than 8) can affect
% the orthogonality

%% Gray levels (discretization levels of the mask)             
discretization = 1; % Variable for the next switch
switch discretization % Gray-level discretized azimuthal angle vector
 case 1 % 1: Evenly-spaced gl phase values.
  gl = 255; % Number of gray levels . Default: 256         
  phaseValues = linspace(0,2*pi,gl); % Discretized phi vector on [0,2*pi]  
  % The sampling interval consists on dividing the range over the gray
  % levels. Similar to the VPL Edgar's discretization formula on the first
  % page of:
  % 1_edgar_2013_High-quality optical vortex-beam generation_E-Rueda_OL.pdf 
  
 case 2 % 2: user-defined gl values
  phaseValues = [0 64 128 255]; % Range of each value: [0,255]
                                % Example: [0 64 128 255]
  % Custom gl vector: the mask will only have these levels
  phaseValues = phaseValues*2*pi/255; % Conversion from gray-levels to 
                                      % phase levels              
  % phase values = gray level values * 2Pi/255. Current range: [0,2pi], but
  % inside f_discretizeMask.m, the values are adjusted to be on [-Pi,Pi] 
end
% Number of gray levels gl = length(phaseValues): calculated inside the
% functions



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% PART 4: MEASUREMENT ADJUSTMENT %%%%%%%%%%%%%%%%%%%%
%% Measurement
% Always saves with the dataformat
% Sweeps all the GL for one tc and then switches to the other ones
% savetype = 1; % 1: as a .mat files or 
                % 2: as dataformat files
dataformat = '.bmp'; % Applies only for savetype = 2
% tcvect = [1 2 3 4 5 6 7 8 9 10]; % Dados por Juan Jose
% glvect = [1 16 24 28 36 56 128 256]; % Dados por Juan Jose
% glvect = [3, 127, 203, 59, 167] % Andres F. Izquierdo: best gl
                                  % with a good system phase response
tcvect = [2]; % Topological charges to be measured
glvect = [2,5]; % Gray level to be measured
% glvect = linspace(2,18,9)
wait = 0; % 10 seconds before measuring as a safety measurement
          % RIGHT NOW 0 FOR DEBUGGING
recordingDelay = 5; % Waits 5 seconds between each mask to be shown
                    % This time is also important so that the camera
                    % bus doesn't overload. Ref: 5
                    % RIGHT NOW 1 FOR DEBUGGING

%% Folder names
if ispc %  Linux (0) or Windows (1)
    pathSep = '\'; % Works for windows 7 and 10
else
    pathSep = '/'; % Works for Linux
end
infoDelim = '-'; % For the data and output folder information
dirDelim= '_'; % "infoDelimiter" is established with the '_' and then
                 % the delimiter for the folder counter should be different
                 % "infoDelim" must always be different than "dirDelim"
analysFldr = 'Analysis'; % Folder name: scripts
dataFlrd = 'Data'; % Folder name: input data  
snapsfldr = 'TestSnapshots'; % Snapshot tests folder (inside dataFlrd)
outFlrd = 'Output'; % Folder name: output data
toolsFldr = 'Tools'; % Folder name: functions
filemanag = 'Files-Folders_Managing'; % Folder with the function 
                                      % f_makeParentFolder, the 1st 
                                      % function that is used in the 
                                      % program
% Folders:
%  Analysis: principal scripts
%  Data: the inputs of the algorithm are the acquired vortex images
%  Data -> DatalogDir: specific measurement folder 
%  Output: processed images or plots
%  Output -> ProcessedDir: specific processed images folder
%  Tools: functions used in the program
%  Tools -> File_managing: where f_makeParentFolder.m and f_addDirectories
%  are



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% PART 5: PROCESSING WITH METRICS %%%%%%%%%%%%%%%%%%%%%%
metricSel = 1; % 
metricProfile = 1; % 1: vertical profile; 2: horizontal profile

%% Diffraction angle units
n = 1; % Air's refractive index
% PP: used with the "camera" variable
M = 1; % No microscope objective is used and there's no magnification
f = 100; % focal length of the lens before the camera [mm]




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% PART 6: ACADEMIC-PURPOSE ASPECTS %%%%%%%%%%%%%%%%%%

% Zernike reconstruction, FT, gradient and a simulation in the free space 
% that is not very depured
%% Optional plots and procedures
FTmask = 0; % Finds the FFT of the mask and plots it: yes(1); no(0)
%%% For FTmask = 1 and abs_ang = 1:
   maskFTlog = 1; % (1)Plots the log10 of the spectrum. (0) normal spectrum                
gradMask = 0; % Finds the gradient of the mask and pltos it: yes(1); no(0)
maskZernReconstr = 0; % Reconstructs the mask with Zernike polynomials and
                      % plots the error  
simBool = 0; % Simulate: yes (1) or no (0)                      

%% Zernike reconstruction parameters
z_ReconstrNumb = 14; % Number of polynomials to use for the reconstruction
% z_pupil: defined in "Zernike polynomials with Noll's convention"               

if simBool == 1
 %% Simulation parameters
 starAmplitude = 1; % ref: 1
 planetAmplitude = 0.8; % ref: 0.1
 pixelSeparation = 0.3; % Linear separation between the bodies; ref: 0.05
 w1 = 0.06; % Beam one width; ref: 0.03
 w2 = w1; % Beam two width; ref: w1/5
 rPupilSize = 0.5; % Input pupil radius: percentage [0,1]; ref: 0.5

 %% Simulation plots to show
 showIin = 1; % Show input intensity
 showPupilin = 1; % Show input pupil
 showFPmag = 1; % Show Fourier plane magnitude
 logscale = 1; % Apply logscale to the Fourier plane magnitude
 showFPphas = 1; % Show Fourier plane phase
 showPhasout = 1; % Show output phase
 showMagout = 1; % Show output magnitude
 showIout = 1; % Show output intensity
end
