function f_ProcessData(measfullpath,refmeasfullpath,ProcessedDir,dataDir,pathSep,infoDelim, ...
dataformat,imgformat,cameraPlane,totalImgs,AiryFactor,metricSel,metricProfile, ...
beepSound,L,NA,PP,mainLyotRadius,measSimulated)

% Inputs:
%   
%
% Outputs:

% Post-processing of the data (application of the metric of the degree of
% extintion)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% INITIALIZATION
%% Processing initialization

%%% Measure the processing time
% Copyright PhD student Jens de Pelsmaeker VUB B-PHOT 2018,Brussels,Belgium
t1_dt = datetime; % store time
disp('Processing started:'); disp(t1_dt)
processedImgname = strcat(ProcessedDir,pathSep,'processed',infoDelim, ...
                          cameraPlane,infoDelim);

if measSimulated == 0 % Real measurement
    %%% Loading all the measurements
    % Explanation: load(directory+filename,variables)
    struct = load(measfullpath); % Loads all the measured images and their info
                            % Two variables are loaded in the structure: 
                            % "expImgs" and "MeasInfo"

    %%% Experimental images
    expMeas = struct.expImgs;

    %%% Experiment information
    measInfo = struct.MeasInfo;

    %%%  Loading the reference measurement
    % The image is read in a uint8 format: integer with values that are
    % normally in [0,255] (8-bit depth or dynamic range)
    refMeas = imread(refmeasfullpath);       
    % since refMeas is a bmp image, it is loaded as uint8

    %%% UINT8 format to Double for the reference image
    % im2double duplicates the precision of the exponent leaving intact the
    % mantisa. It as floating-point format that normalizes the images and this
    % operation is made on each RGB channel. rgb2gray does a similar operation
    % but scaling to a gray scale, where it  converts RGB values to grayscale values
    % by forming a weighted sum of the R, G, and B components:
    % 0.2989 * R + 0.5870 * G + 0.1140 * B 
    refMeas = im2double(refMeas);
    
   
    
    % Lyotimg = refMeas;  % TO BE USED FOR LYOT METRICS                                                                                                   
    
else % Simulated measurement
   %% Example images to process
%     Lyotimg = imread(strcat(dataDir,pathSep,'0_ExampleData',pathSep,'data_ref_1.bmp')); % Lyot image
%     Lyotimg = rgb2gray(Lyotimg);  
    refMeas = imread(strcat(dataDir,pathSep,'0_ExampleData',pathSep,'data_ref_2.bmp')); % PSF reference
    expMeas = {0,0}; % Cell initialization
    expMeas{1} = imread(strcat(dataDir,pathSep,'0_ExampleData',pathSep,'data_ref_3.png')); % PSF measurement 1
    expMeas{2} = imread(strcat(dataDir,pathSep,'0_ExampleData',pathSep,'data_ref_4.png')); % PSF measurement 2
    measInfo = {'data_ref_3','data_ref_3'};                    
    totalImgs = 2; % For the case of the simulated measurement
    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TO UNIFY WITH DEFINE SPACE
%% Find the center of the Lyot image
% This was already done in f_DefineSpace.m
% PP = PP*1e-6; % um to m

%%% Lyot's spot size (main radius)
%%%%% Lyot Intensity Feedback coordinates
% drawing = false;

% [~,mainLyotRadius,~] = f_findCircleShapedIntensity(Lyotimg,drawing);
% mainLyotRadius = round(mainLyotRadius); % Pixels

%%%%% System Pixel Size:
% not used since the spot size is used intead for the falco lamda over D (physical scaling)
% PP = PP*1e-6; % um to m
% lensDiameter = 2*apRad*1e-2; % cm to m
% [apRadpix] = f_computePupilPixelSize(mainLyotRadius,PP,lensDiameter);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TO UNIFY WITH DEFINE SPACE



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Coordinates: Lambda/D, pixels and arcseconds
%% Find the center of the PSF image (with a binarization)                                                                   
[~,~,aproxCenter,aproxRadius] = f_approximateSpotSize(refMeas);
 
% MAYBE: Find center of the PSF image (peaks)  [for more robustness]

%% Read the approximate center of the PSF reference
midX = aproxCenter(2);
midY = aproxCenter(1);

%%% OLD: center of the image but not the spot's center
% midX = round((maxX+1)/2); % x mid point
% midY = round((maxY+1)/2); % y mid point
                        
%% Cartesian coordinates with pixel units
% ORIGINAL
[ySize, xSize] = size(refMeas); % All images assumed of the same size as the refmeas
halfX = xSize/2;
halfY = ySize/2;

% Old: 
% halfX = floor((ySize+1)/2); % x mid point
% halfY = floor((xSize+1)/2); % y mid point

% To delete:
% [ySize,xSize] = size(Lyotimg);
% halfX = xSize/2;
% halfY = ySize/2;

%%  lambda/D factor falco-matlab reference
% It is scalled with respect to the jinc zeros

NpadX = xSize; % Camera's x pixel size
NpadY = ySize; % Camera's y pixel size

xlamOverD = NpadX/(2*mainLyotRadius);
ylamOverD = NpadY/(2*mainLyotRadius);

% Centering shifting to the spot location
centerShiftX = (midX-halfX); 
centerShiftY = (midY-halfY);

% Coordinates' origin set to the spot's center
xpixcenterd= (-halfX:halfX-1) - (centerShiftX - 1); % The center has to be shifted 1
ypixcenterd = (-halfY:halfY-1) - (centerShiftY - 1); % The center has to be shifted 1

% xangL_Dfalco = xpixcenterd/xlamOverD; % Astronomer's physical scaling of pixels
% yangL_Dfalco = ypixcenterd/ylamOverD; % Astronomer's physical scaling of pixels

%% Lambda over D scaling with the experimental spot size
% Pixel's size is scalled to the spot's size
xangL_Dexpairy = xpixcenterd/(aproxRadius); % aproxRadius/2 makes it the diameter
yangL_Dexpairy = ypixcenterd/(aproxRadius);

%% Lambda over D scaling with the experimental spot size
% Pixel's size is scalled to the first Bessel's center
% AiryFirstZero = 1.22; % First zero of the cylindrical Bessel function of 
%                         % first kind and zeroth order
% xangL_Dexpairyzero = xangL_Dexpairy*AiryFirstZero;
% yangL_Dexpairyzero = yangL_Dexpairy*AiryFirstZero;

%% Symmetric pixels (sign type)
% xpixsym = -halfX : halfX - 1; % pixels
% ypixsym= -halfY :  halfY - 1; % pixels      

%% Unitary pixels (scaled heavyside)
% xpix = 1:xSize; % Pixels start in 1
% ypix = 1:ySize; % Pixels start in 1

%% Cartesian coordinates with the lambda/D scaling (diffraction angle)
                                                                                                                             %  AiryFactor is an input
% xangL_Dtheoric = f_scalePix2DiffAng(xpixcenterd,AiryFactor);
% yangL_Dtheoric = f_scalePix2DiffAng(ypixcenterd,AiryFactor);

%% Definitive Lambda over D vector
% xangL_D
% yangL_D

%% Cartesian coordinates with the arcsecond scaling (diffraction angle)
% angArcs = f_LambdaDToarcsec(xangL_D);
% yangArcs = f_LambdaDToarcsec(yangL_D);

%% Cartesian coordinates selector
x = xangL_Dexpairy;
y = yangL_Dexpairy;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Profiles for the metrics
%% Reference Profile
%%% Metric-specific default parameters for the profile
oneSideProfile = 1; % Specifically needed for this metric. Ref: 1
shiftCart = [0,0]; % midX,midY already account for the shift

%%% Find the reference profile
[xrefProx,yrefProf,HprofRef,VprofRef,~,~] = f_makeImageProfile(x,y,midX,midY, ...
           refMeas, shiftCart, oneSideProfile);
[averRefProfile] = f_getAverageRadialProfile(refMeas,[ySize, xSize],aproxCenter);

%% Profile of the measurements
Hprofmeas = cell(1,totalImgs);
Vprofmeas = cell(1,totalImgs);
averMeasProfile = cell(1,totalImgs); 

for idxgral = 1:totalImgs
    [~,~,Hprofmeas{idxgral},Vprofmeas{idxgral},~,~] = f_makeImageProfile(x,y,midX,midY, ...
       expMeas{idxgral},shiftCart, oneSideProfile);
    [averMeasProfile{idxgral}] = f_getAverageRadialProfile(expMeas{idxgral},[ySize, xSize],aproxCenter);
end

%% Measurement profile choice
switch metricProfile
  case 1 % Vertical profile
    radialIntensityMeas = Vprofmeas; % One-sided
    
  case 2 % Horizontal profile
    radialIntensityMeas = Hprofmeas; % One-sided
    
  case 3 % Radial averaged profile
    radialIntensityMeas = averMeasProfile;
    
  otherwise
    error('"metricProfile" must be either 1 or 2');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Metric application
%% Processsing of the profiles
for idxgral = 1:totalImgs
    [measEEFcurves{idxgral}] = f_calculateEEF(radialIntensityMeas{idxgral},cartcoord,titprof,tit,xlab);
end

for idxgral = 1:totalImgs                                                                                                                    % ACTUALLY, CALCULATE ALL THE METRICS BUT ONLY PLOT THE WANTED ONES
  switch metricSel
    case 1 
     %% Throughput: Encircled Energy Factor metric
      % Camera: PSF.
      tit = 'Encircled Energy Distribution of Intensity';
      [energy] = f_calculateEEF(radialIntensityRef);
     
     %% Throughput gradient
                                                                                                                                                                            % DO A FUNCTION ONLY IF THIS WILL BE USEFULL
%       tit = 'Throughput gradient';                                               
%       normIntensity = radialIntensityMeas./max(radialIntensityMeas);
%       GradEnergy = gradient(normIntensity); % gradient [returns n elements] or diff [returns n-1 elements]                            % OR GRAD energy ?
      
     %% Plot of the gradient of the EEF and its corresponding intensity pattern                                                              
      %       yyaxis left (OLD method)
      %       h1 = plot(cartcoord,GradEnergy);
      %       yyaixs right
      %       h2 = plot(cartcoord,normIntensity);
      figure('color','white');
      [hAxes,hLine1,hLine2] = plotyy(cartcoord,GradEnergy,cartcoord,normIntensity);
      set(hLine1, 'Color','b','LineStyle','-'); set(hLine2, 'Color','r','LineStyle','-');
      set(hAxes(1),'YColor','b'); set(hAxes(2),'YColor','r'); set(hAxes,'FontSize',12,'FontWeight','bold');
      axis(hAxes(1), 'square'); axis(hAxes(2), 'square'); box(hAxes(1), 'on');
      title(strcat(tit,{' '},titprof)); grid on;
      xlabel('Angular position [\lambda/D]'); 
      ylabel(hAxes(1),'Gradient of the relative throughput','FontSize',12,'FontWeight','bold');
      ylabel(hAxes(2),'Intensity pattern','FontSize',14,'FontWeight','bold'); grid on;
      lgdHandler = legend({'Gradient of the Encircled Energy Factor', 'Intensity'}); lgdHandler.FontSize = 10;
      
     %% Power suppresion in the airy disk
      tit ='Power suppresion in the airy disk';    
      % Needs case 1 as input, BUT FOR MULTIPLE GLs to do one cycle. A plot
      % is done per TC !
      
    case 2 % SNR
      tit = 'Signal-to-Noise Ratio';
      tit1 = strcat(tit,{' '},titprof);
      tit2 = strcat('Raw contrast',{' '},titprof);
      [SNR] = f_calculateSNR(radialIntensityMeas,radialIntensityRef, ...
      cartcoord,tit1,tit2,xlab);
    case 3 % MSE
      tit = 'Mean Squared Error';
        
    otherwise
      error('Select a valid metric');
  end
  

  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%  Save and finish the processing

  %% Saving                                                                      
  %                                                                                                                                             FUTURE PLOT SAVING
  processedImgfullpath = strcat(processedImgname,measInfo{idxgral});
  % Explanation: saveas(variable,directory+filename,extension)
  saveas(gcf,strcat(processedImgfullpath),imgformat); % Saves the last shown figure
  
  % OLD:
  %   % Explanation: imwrite(variables,directory+filename+extension)
  %   imwrite(expMeas{idxgral}, strcat(processedImgfullpath,dataformat));
  
end

%% End of the processing
% Author: PhD student Jens de Pelsmaeker VUB B-PHOT 2018, Brussels, Belgium
% MATLAB built in:
t2_dt = datetime;
disp('Processing finished:'); disp(t2_dt)
time = t2_dt - t1_dt;
disp('Processing took: '); % datestr(time,'SS') ' seconds'])
disp(time);

%% End notification
N = 4; % Number of beeps for the processing
f_EndBeeps(N,beepSound);

%% Ask to leave figures open
answer = questdlg('Do you want to close all the figures?','Processing finished','yes','no','no');
if strcmp(answer,'yes') % Compare string
     close all;         
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Airy size estimatives
%% Theoretical: diffraction-limited systems
% NA: numerical aperture of the lens
% L: wavelength in um
PP = PP*1e-6; % um to m % Measurement camera's PP

AiryDiskSpatial = 0.61*L/NA; % um
AiryDiskSpatialX = AiryDiskSpatial;
AiryDiskSpatialY = AiryDiskSpatial;
disp(strcat('Theoretical Airy"s radius in um:', {' '}, num2str(AiryDiskSpatialX), ' (x) and', {' '},  num2str(AiryDiskSpatialY), ' (y)'));

AiryDiskPixX = round(AiryDiskSpatialX/PP); % PP: camera's pixel pitch
AiryDiskPixY = round(AiryDiskSpatialY/PP);
disp(strcat('Theoretical Airy"s radius in pix:', {' '},  num2str(AiryDiskPixX), ' (x) and', {' '},  num2str(AiryDiskPixY), ' (y)'));

% In reality, there are aberrations and the real radius can be of
% about 60 times the theoretical one

%% Airy radius measured from the reference image
AiryDiskPixX = aproxRadius; % Just an example
AiryDiskPixY= aproxRadius; % Just an example
AiryDiskSpatialX = AiryDiskPixX*PP; % PP: camera's pixel pitch
AiryDiskSpatialY = AiryDiskPixY*PP;
disp(strcat('Estimated Airy"s radius in um:', {' '},  num2str(AiryDiskSpatialX), ' (x) and', {' '},  num2str(AiryDiskSpatialY), ' (y)'));
disp(strcat('Estimated Airy"s radius in pix:', {' '},  num2str(AiryDiskPixX), ' (x) and', {' '},  num2str(AiryDiskPixY), ' (y)'));

%% Airy radius from the EEC factor
% When it is the 70%

%% Pixel airy radius
% AiryMultiplicityX = xSize/AiryDiskPixX; % Number of airy disks 
% AiryMultiplicityY = ySize/AiryDiskPixY;

%% Lambda over D scaled coordinates OLD METHOD
% xangL_D = xpix/(AiryMultiplicityX*8);
% yangL_D =  ypix/(AiryMultiplicityY*8);

end