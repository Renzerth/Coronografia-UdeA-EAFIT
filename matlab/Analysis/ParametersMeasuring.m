% COMMENTS ARE ONLY MADE IN Parameters.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% PART 1: GENERAL ADJUSTMENTS %%%%%%%%%%%%%%%%%%%%%%%
%% Algorithm sections
meas = 1; 
    measSimulated = 0; 
    measDebug = 0;
    beepSound = 1; 
proc = 0; 
    loadMeas = 1; 
    measFoldName = '04-Apr-2019-No-SLM-VPL-mask-tcs-1-gls-1';

%% General algorithm parameters: coordinates, plots, screens and mask type
precision = 3; 
abs_ang = 2;
maskSel = 0; 
plotMask = 2;        
                     
                     
              

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%              
%%%%%%%%%%%%%%%%%%%%%%% PART 2: HARDWARE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
slm = 'Pluto'; 
switch slm
  case 'Pluto'
    %% SLM parameters (reflection)
    sSupport = min([1.536 0.864]); 
    maxNumPix = max([1920 1080]); 
    pixSize = 8; 
    scrnIdx = 3; 
    
  case 'LC2002'
    %% SLM parameters (transmision)
    sSupport = min([2.66 2.00]); 
    maxNumPix = max([800 600]); 
    pixSize = 32; 
    scrnIdx = 2; 
    
  case 'No-SLM'
   sSupport = 1; 
   maxNumPix = max([1920 1080]); 
   pixSize = 8; 
   scrnIdx = 1; 
   
  otherwise
    warning('Please select a valid SLM');
end 

%% SLM positionining calibration, coordinates and type of truncation
MaskPupil = 1;
coordType = 2; 
    k = 9; 
    sSize = 2^k-1;  
    MaxMask = 1; 
    circularMask = 1; 
    shiftMask = 2; 
     switch shiftMask
         case 1 
              shiftCart = [-10, -20]; 
                        
         case 2 
             shiftCart = [0, 0];
     end
              
%% Measurement camera selection and parameters
camera = 'DMK42BUC03'; 

switch camera  
  case 'DMK42BUC03' 
    cameraID = 2;
    exposure = 1/250; 
    fps = [];
    format = 'Y800 (1280x960)';
    cameraPlane =  'Lyot';
    PP = 3.75;
 
  case 'DMK23U445' 
    cameraID = 1;
    exposure = 1/300; 
    fps = '15.00'; 
    format = 'Y800 (1280x960)'; 
    cameraPlane = 'PSF';
    PP = 3.75;
    
  case 'DMK41BU02.H' 
    cameraID = 3;  
    exposure = 1/1e3; 
    fps = '15.0000';
    format = 'Y800 (1280x960)';    
    cameraPlane = 'notusedhere';
    PP = 4.65; 
end

%% Image capture for the measurement debug
filename = 'test'; 
imgformat = 'png'; 
previewBool = 1;
loghist = 0; 
                    
                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          
%%%%%%%%%%%%%%%%%%%%% PART 3: PHASE MASKS PARAMETERS %%%%%%%%%%%%%%%%%%%%%%
%% Parameters: Laguerre-Gauss, spiral phase mask and general masks
L = 0.6328; 
tc = 2; 
s = 1;
ph0 = 0;
    binMask = 0;
    binv = 0; 
    normMag = 1; 
          
%% Parameters: Laguerre-Gauss
p = 4; 
WsizeRatio = 20; 
         
%% Parameters: VPL Phase mask, 
f_FR = 1.3; 
AreaSLM = maxNumPix*pixSize^2; 
minf_FRum = AreaSLM/L; 
scaleFactor = 1e-6;
minf_FR = minf_FRum*scaleFactor;
if f_FR < minf_FR 
  error(['VPL criterium not fulfilled, establish f_FR bigger than ' ...
         num2str(minf_FR)]);
end

%% Parameters: Elliptic Gaussian Vortex
bcst = 0.3;

%% Parameters: Fork Phase
frkTyp = 1;
period = 0.1; 
T0 = 1; 
Aalpha = pi; 
Angalp = pi/2; 
Angbet = 0; 

%% Parameters: Zernike polynomials with Noll's convention
z_coeff = [4 5]; 
z_a = 2.5; 
z_pupil = 1; 
z_disp_wrap = 1; 
z_plot = 0; 

%% Gray levels (discretization levels of the mask)         
discretization = 1; 
switch discretization 
 case 1 
  gl = 256; 
  phaseValues = linspace(0,2*pi,gl);
  
 case 2 
 phaseValues = [0 255];
 
 phaseValues = phaseValues*2*pi/255; 
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% PART 4: MEASUREMENT ADJUSTMENT %%%%%%%%%%%%%%%%%%%%
%% Measurement
dataformat = '.bmp'; 
tcvect = [1 2 0]; 
glvect = [255]; 
waitbeforemeas = 2; 
recordingDelay = 0; 

%% Reference measurement     
whiteblackref = 0; 

%% Folder names
if ispc 
    pathSep = '\'; 
else
    pathSep = '/'; 
end
infoDelim = '-'; 
dirDelim= '_';          
analysFldr = 'Analysis'; 
dataFlrd = 'Data';
snapsfldr = 'TestSnapshots'; 
outFlrd = 'Output';
toolsFldr = 'Tools';
lastmeas = 'LastMeasName'; 
filemanag = 'Files-Folders_Managing'; 
      
                                      


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% PART 5: PROCESSING WITH METRICS  AND SELF CENTERING ALGORITHM %%%%%%%%

metricSel = 1;
metricProfile = 1; 

%% Optical system parameters
n = 1; 
M = 10;
f = 100; 
scaleFactor = 1e3; 
f = f*scaleFactor; 
AiryFactor = n*PP*M/f; 
NA = 0.1;
apRad = 2.54; 

                       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% PART 6: ACADEMIC-PURPOSE ASPECTS %%%%%%%%%%%%%%%%%%

%% Optional plots and procedures
FTmask = 0; 
abs_angFT = 1; 
   maskFTlog = 1;            
gradMask = 0; 
maskZernReconstr = 0; 
simBool = 0;               

%% Zernike reconstruction parameters
z_ReconstrNumb = 14;       

if simBool == 1
 %% Simulation parameters
 starAmplitude = 1;
 planetAmplitude = 0.8;
 pixelSeparation = 0.3; 
 w1 = 0.06; 
 w2 = w1; 
 rPupilSize = 0.5; 

 %% Simulation plots to show
 showIin = 1;
 showPupilin = 1; 
 showFPmag = 1; 
 logscale = 1;
 showFPphas = 1;
 showPhasout = 1;
 showMagout = 1; 
 showIout = 1; 
end