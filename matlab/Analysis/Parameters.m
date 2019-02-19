%%%%%%%%%%%%%%%%%%%%%%% PART 1: GENERAL ADJUSTMENTS
%% Algorithm sections
measDebug = 0; % Debugging before actually measuring. Displays the default 
               % phase mask and shots a photo with the camera
meas = 0; % Measure: yes (1) or no (0)
beepSound = 0; % Beep sound when measurement finishes. Only works when 
               % meas = 1
slmselect = 1; % 1: Pluto (reflection); 2: LC2002 (transmission)
sim = 0; % Simulate: yes (1) or no (0)

%% General algorithm parameters
coordType = 2; % 
  k = 10; % Bits for grey levels; 2^k is the resolution (size of x and y)
        % Default: 10        
precision = 3; % Precision of displayed results: significative digits (3)
showM = 1; % Plot the individual mask inside "PhaseMaskSel.m": no(0)-yes(1)
           % analog to "plotMask" on the SLM Position section
maskSel = 0; % Phase mask selection:
% 0: Helicoidal mask: SPP or DSPP depending on gl
% 1: Laguerre-Gauss beams: amplitude or phase
% 2: VPL: Vortex Producing Lens = Helicoidal + Fresnel lens
% 3: Elliptic Gaussian beam phase mask
% 4: Fork phase masks
% ---- NOT USED:
% 5: Zernike (aberrations)
% 6: Laguerre-Gauss + Zernike
% 7: Hermite-Gauss beams NOT DONE
% 8: Mutliple vortices NOT DONE
% 9: Sum of spiral phase masks NOT DONE
% 10: Gerchberg-Saxton NOT DONE
% otherwise: Unitary




%%%%%%%%%%%%%%%%%%%%%%% PART 2: HARDWARE
% screenIndex: screen number selector. In [0,N-1] with N the # of screens
if slmselect  == 1
    %% SLM parameters (reflection)
    % spaceSupport = 1; % Unitary space: spaceVector = -1:2/(Ssize-1):1;
    spaceSupport = min([0.864 1.536]); % Size of the SLM window in cm:
                                         % 1.536cm x 0.864cm
    maxNumPix = max([1920 1080]); % Maximum number of pixels on the SLM 
                                  % (either horizontal or vertical); SLM's 
                                  % resolution in pixels: 1920 x 1080 
    pixSize = 8; % SLM pixel's size in um
    scrnIdx = 2; % Screen number selector
else
    %% SLM parameters (transmision)
    spaceSupport = min([2.66 2.00]); % Same as the reflection SLM
    maxNumPix = max([800 600]); % Same as the reflection SLM
    pixSize = 32; % Same as the reflection SLM in um
    scrnIdx = 2; % Screen number selector
end 

%% SLM positionining calibration
shiftBool = 0; % Shift activated (1)[SLM displaying] or deactivated (0)
               % [exporting masks]. shiftCart = [yshift,xshift]
shiftCart = [0,0]; % Percentages of movement of the total size of the
                   % mask (cartesian coordinates convention)
                   % Calibrated with: s = +1; ph0 = 0, tc = 1; 
plotMask = 2; % Allows to plot the final mask, as it can be a combination 
              % of the previous ones
              % 0: no plot;
              % 1: on the screen
              % 2: on the SLM
              % 3: on the screen but surface-plot type
              
%% Camera selection
camera = 'DMK23U445';
% 'DMK23U445' : PSF plane. Exposure: [,]
% 'DMK42BUC03' : Lyot plane. Exposure: [1/1e4,1]
% 'DMK41BU02.H' : not used here.
exposure = 1/1e3; % Analog parameter. 
% For 'DMK42BUC03':
% -Default: 0.0556 or 1.282
% -PH50um: 0.0030 -- Y0.0256 || PH25um: 0.0083 -- Y0.0556 ||
% -PH15um: --Y0.1429 || 0.0227 || PH10um: 0.0435 -- Y0.200 || 
% -PH5um: 0.363 -- Y3.099
format = 'Y800 (1280x960)'; % Best: 'Y800 (1280x960)'
                            % Another option: 'RGB24 (1024x768)'

%% Image capture
filename = 'test'; % Name of the capture one wants to take
imgformat = '.png'; % Format with period. mat, bmp, png, jpg



%%%/%%%%%%%%%%%%%%%%% PART 3: DEFAULT PHASE MASKS
%% Parameters: Laguerre-Gauss, spiral phase mask and general masks
L = 0.6328; % Laser wavelength [um]. Used in Zernike and VPL masks
abs_ang = 2; % Magnitude (1) or phase (2) plot
tc = 2; % Topological charge (integer bigger or equal to one)
        % tc = Azimuthal index m for LG. Fractional tc result on phase
        % patterns of Hermite-Gauss (maybe just a coincidence)
s = 1; % Sign of mask (+1 or -1); reverses the imprinted OAM 
ph0 = 0; % Initial phase of the angle [radians]; reference +pi from
         % normal zero of trig circle and same rotation convention.
         % This corresponds to a normal rotation of the mask for stethic
         % reasons and shouldn't affect the results. Only affects if the
         % vortex is no fully centered
binMask = 0; % Binarizes the mask w.r.t the max/min of the phase (boolean)
       
%% Gray levels (discretization levels of the mask)
% Dynamic range = maxGrayDepth - minGrayDepth
mingl = 0; % Minimum gray level depth. Ref: 0
maxgl = 255; % Maximum gray level depth. Ref: 255
levShft = 0; % Ref: 0. Seems to be non-linear or better not to use it
             % Corresponds to the brightness or constant shift of the gl's
gl = 2; % Number of grey levels (normally 256). Must be smaller than 5
          % the dynamic range = maxGrayDepth-minGrayDepth

%% Parameters: Laguerre-Gauss
p = 5; % Number of radial nodes. If p=0, normal helicoid masks are obtained
       % If they are used and tc=0(m=0); binary masks are obtained
       % Even p; rings are ones. Odd p; rings are zeroes. Use mask = mask'
W = 100; % Width of the modes; for LG; ref: 100
binv = 1; % Binary inversion of the mask: yes(1); no(0). It is only applied 
          % when tc is zero. It is usefull to be applied for odd p. 
norm = 0; % Normalize magnitude and phase (to unity). yes(1); no(0)         

%% Parameters: VPL Phase mask, 
% f_FR: Fresnel lens focal distance or diffractive lens phase focal length
f_FR = maxNumPix*pixSize^2/L; % Criterium to determine the MINIMUM f_FR
% From: 2_edgar_2015_Generation_Optical_Vortices_Binary_Vortex_Lenses.pdf


%% Parameters: Elliptic Gaussian Vortex
bcst = 0.1; % Ellipticity. cy/cx = 1/alpha. Ref: .1, .2, .4, .6, .8 and 1
            % Out of theory: when beta >> 1, the mask tends to be binary
            % When beta ~ 0, the mask tends to be trinary

%% Parameters: Fork Phase
frkTyp = 2; % 1: smooth transition (amplitude); 2: phase jump transition (phase)
            % 1 = 2 when binarized
period = 0.1; % Period of the grating (fringe spacing). Ref: 0.1
              % Twice its inverse is the number of line dislocations
              % So the frequency is f = 2/period
              
%%% Smooth transition
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

%%% Phase jump transition
Angalp = 0; % Diffraction angle of horizontal direction (x) [radians]
Angbet = pi/2; % Diffraction angle of vertical direction (y) [radians]
% Here, more than pi/2 seems not to work very well
% Range: [-pi/2,pi/2]
% The "*" in the 'Smooth transition' case also apply here




%%%%%%%%%%%%%%%%%%%%%%% PART 4: MEASUREMENT ADJUSTMENT
%% Measurement
% tcvect = [1 2 3 4 5 6 7 8 9 10]; % Dados por Juan Jose
% glvect = [1 16 24 28 36 56 128 256]; % Dados por Juan Jose
% glvect = [3, 127, 203, 59, 167] % Andres F. Izquierdo: best gl
                                  % with a good system phase response
tcvect = [1]; % Topological charges to be measured
glvect = [255]; % Gray level to be measured
wait = 0; % 10 seconds before measuring as a safety measurement
recordingDelay = 1; % Waits 5 seconds between each mask to be shown
                    % This time is also important so that the camera
                    % bus doesn't overload. Ref: 5
                    % RIGHT NOW 1 FOR DEBUGGING

%% Folder names
analysFldr = 'Analysis'; % Folder name: scripts
dataFlrd = 'Data'; % Folder name: input data           
outFlrd = 'Output'; % Folder name: output data
toolsFldr = 'Tools'; % Folder name: functions



              
%%%%%%%%%%%%%%%%%%%%%%% PART 5: ACADEMIC PURPOSES: Zernike, simulation
%% Parameters: Zernike
% L and gl are used here
z_coeff = [1, 0.1]; % Zernike coeffient vector (see f_Zernike_Mask.m)
a = 60; % Arbitrary constant; the bigger, the more intense; ref: a=20
frac = 0.125; % To adjust the wrapped phase; ref: 0.125
pupil = 1; % Pupil relative size: [0,1]; like a percentage
disp_wrap = 0; % (0): Original; (1): wrapped mask on [-pi,pi] 
plot_z = 0; % plot with Zernike builder: yes(1); no(0)

%% Optional plots and procedures
FTmask = 0; % Finds the FFT of the mask and plots it: yes(1); no(0)
gradMask = 0; % Finds the gradient of the mask and pltos it: yes(1); no(0)
maskZernReconstr = 0; % Reconstructs the mask with Zernike polynomials and
                      % plots the error
                      

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

