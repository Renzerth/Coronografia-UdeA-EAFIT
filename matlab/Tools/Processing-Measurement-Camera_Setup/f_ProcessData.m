function f_ProcessData(measfullpath,refmeasfullpath,ProcessedDir,dataDir,pathSep,infoDelim, ...
dataformat,imgformat,cameraPlane,totalImgs,AiryFactor,metricSel,metricProfile, ...
shiftCart,beepSound,L,NA,PP,apRad)

% Inputs:
%   
%
% Outputs:

%% Post-processing of the data (application of the metric of the degree of
%%% extintion)

%% Measurements initializing
% Copyright PhD student Jens de Pelsmaeker VUB B-PHOT 2018,Brussels,Belgium
t1_dt = datetime; % store time
disp('Processing started:'); disp(t1_dt)
processedImgname = strcat(ProcessedDir,pathSep,'processed',infoDelim, ...
                          cameraPlane,infoDelim);

%% Loading all the measurements
% Explanation: load(directory+filename,variables)
struct = load(measfullpath); % Loads all the measured images and their info
                        % Two variables are loaded in the structure: 
                        % "expImgs" and "MeasInfo"
                        
%%  Loading the reference measurement
% The image is read in a uint8 format: integer with values that are
% normally in [0,255] (8-bit depth or dynamic range)
refmeas = imread(refmeasfullpath);       

%% UINT8 format to Double for the image
% im2double duplicates the precision of the exponent leaving intact the
% mantisa. It as floating-point format that normalizes the images and this
% operation is made on each RGB channel. rgb2gray does a similar operation
% but scalling to a gray scale, where it  converts RGB values to grayscale values
% by forming a weighted sum of the R, G, and B components:
% 0.2989 * R + 0.5870 * G + 0.1140 * B 
refmeas = im2double(refmeas);

%% Find center of the Lyot image
% This was already done in f_DefineSpace.m

%% Find center of the PSF image (binarization)
[xangL_D,yangL_D,regionCentroid,aproxRadius] = f_approximateSpotSize(refmeas);

%% Find center of the PSF image (peaks)                                                  MISSING!!
                        
%% Cartesian coordinates with pixel units
% [ySize, xSize] = size(refmeas); % All images assumed of the same size as the refmeas
% xpix = 1:xSize; % Pixels start in 1
% ypix = 1:ySize; % Pixels start in 1

%% Theoretical: diffraction-limited systems
% NA: numerical aperture of the lens
% L: wavelength in um
AiryDiskSpatial = 0.61*L/NA; % um
AiryDiskSpatialX = AiryDiskSpatial;
AiryDiskSpatialY = AiryDiskSpatial;
disp(strcat('Theoretical Airy"s radius in um: ', num2str(AiryDiskSpatialX), ' (x) and ', num2str(AiryDiskSpatialY), ' (y)'));

AiryDiskPixX = AiryDiskSpatialX/PP; % PP: camera's pixel pitch
AiryDiskPixY = AiryDiskSpatialY/PP;
disp(strcat('Theoretical Airy"s radius in um: ', num2str(AiryDiskPixX), ' (x) and ', num2str(AiryDiskPixY), ' (y)'));

% In reality, there are aberrations and the real radius can be of
% about 60 times the theoretical one

%% Airy radius measured from the reference image
AiryDiskPixX = aproxRadius; % Just an example
AiryDiskPixY= aproxRadius; % Just an example
AiryDiskSpatialX = AiryDiskPixX*PP; % PP: camera's pixel pitch
AiryDiskSpatialY = AiryDiskPixY*PP;
disp(strcat('Estimated Airy"s radius in um: ', num2str(AiryDiskSpatialX), ' (x) and ', num2str(AiryDiskSpatialY), ' (y)'));
disp(strcat('Estimated Airy"s radius in um: ', num2str(AiryDiskPixX), ' (x) and ', num2str(AiryDiskPixY), ' (y)'));

%% Airy radius from the EEC factor
% When it is the 70%

%% Pixel airy radius
% AiryMultiplicityX = xSize/AiryDiskPixX; % Number of airy disks 
% AiryMultiplicityY = ySize/AiryDiskPixY;

%% Symmetric pixels

% xpix = -xSize/2 : 1 : xSize/2 - 1; % pixels
% ypix= -ySize/2 : 1 : ySize/2 - 1; % pixels                   % MAYBE USE MIDX,MIDY

% xangL_D = xpix/(AiryMultiplicityX*8);
% yangL_D =  ypix/(AiryMultiplicityY*8);

%% Cartesian coordinates with the lambda/D scaling (diffraction angle)
% xangL_D = f_scalePix2DiffAng(xpix,AiryFactor);  
% yangL_D = f_scalePix2DiffAng(ypix,AiryFactor);

%% Cartesian coordinates with the arcsecond scaling (diffraction angle)
% xangArcs = f_LambdaDToarcsec(xangL_D);
% yangArcs = f_LambdaDToarcsec(yangL_D);

%% Lyot's spot size (main radius)
PSFimg = imread(strcat(dataDir,pathSep,'data_ref_2.bmp'));
Lyotimg = imread(strcat(dataDir,pathSep,'data_ref_1.bmp'));
Lyotimg = rgb2gray(Lyotimg); 

drawing = false;
PP = PP*1e-6; % um to m
apRad = apRad*1e-2; % cm to m
[~,mainLyotRadius,~] = f_findCircleShapedIntensity(Lyotimg,drawing);
mainLyotRadius = round(mainLyotRadius);
[apRadpix] = f_computePupilPixelSize(mainLyotRadius, PP,apRad);

%% Find center of the PSF image (binarization)
[xangL_D,yangL_D,regionCentroid,aproxRadius] = f_approximateSpotSize(PSFimg);

midX = regionCentroid(2);
midY = regionCentroid(1);
% midX = round((maxX+1)/2); % x mid point
% midY = round((maxY+1)/2); % y mid point

%%  lambda/D factor falco-matlab reference
% It is scalled with respect to the jinc zeros

[ySize,xSize] = size(Lyotimg);
halfX = xSize/2;
halfY = ySize/2;

NpadX = xSize; % Camera's x pixel size
NpadY = ySize; % Camera's y pixel size

xlamOverD = NpadX/(2*mainLyotRadius);
ylamOverD = NpadY/(2*mainLyotRadius);

centerShiftX = (midX-halfX);
centerShiftY = (midY-halfY);

xvals_fp = (-halfX:halfX-1) - centerShiftX;
yvals_fp = (-halfY:halfY-1) - centerShiftY;

xpix = xvals_fp/xlamOverD;
ypix = yvals_fp/ylamOverD;

%% Metric-specific default parameters
xlab = 'Angular position [\lambda/D]';
ylab = 'Angular position [\lambda/D]';
plotData = 1; % Shows the profile lines. Ref: 1
plotH = 0; % Not needed for the metric. Ref: 0
plotV = 0; % Not needed for the metric. Ref: 0
oneSideProfile = 1; % Specifically needed for this metric. Ref: 1
tol = 0; % 0: no need to symmetrically truncate the profile. Ref: 0
shiftCart = [0,0]; % midX,midY already accounts for the shift
tit = 'tit';

%% Reference Profile
[x,y,HprofRef,VprofRef,~,~] = f_makeImageProfile(xpix,ypix,midX,midY, ...
           PSFimg,tol,shiftCart,tit,xlab,ylab,plotData,plotH,plotV,oneSideProfile);

%% Reference Profile choosing 
switch metricProfile
    case 1 % Vertical profile
        radialIntensityRef = VprofRef; % One-sided
        cartcoord = y;
        titprof = '(vertical profile)';
        
    case 2 % Horizontal profile
        radialIntensityRef = HprofRef; % One-sided
        cartcoord =x;
        titprof = '(horizontal profile)';
        
    otherwise
        error('"metricProfile" must be either 1 or 2');
end

for idxgral = 1:totalImgs
%% Profile of the measurements
[xpix,ypix,Hprofmeas,Vprofmeas,~,~] = f_makeImageProfile(xpix,ypix,midX,midY, ...
   struct.expImgs{idxgral},tol,shiftCart,tit,xlab,ylab,plotData,plotH,plotV,oneSideProfile);
       
%% Profile choosing measurements
switch metricProfile
    case 1 % Vertical profile
        radialIntensityMeas = Vprofmeas; % One-sided
        
    case 2 % Horizontal profile
        radialIntensityMeas = Hprofmeas; % One-sided
        
    otherwise
        error('"metricProfile" must be either 1 or 2');
end

  %% Processsing of the image
                                                                                                                      % ACTUALLY, CALCULATE ALL THE METRICS BUT ONLY PLOT THE WANTED ONES
  switch metricSel
    case 1 
     %% Throughput: Encircled Energy Factor metric
      % Camera: PSF.
      tit = 'Encircled Energy Distribution of Intensity';
      [energy] = f_calculateEEF(radialIntensityRef,cartcoord, tit,xlab);
     
     %% Throughput gradient
                                                                                                                                                                            % DO A FUNCTION ONLY IF THIS WILL BE USEFULL
      tit = 'Throughput gradient';                                               
      normIntensity = radialIntensityMeas./max(radialIntensityMeas);
      GradEnergy = gradient(normIntensity); % gradient [returns n elements] or diff [returns n-1 elements]                            % OR GRAD energy ?
      
     %% Plot of the gradient of the EEF and its corresponding intensity pattern                                                                                                                                                    %
                                                                                                                                                 
%       yyaxis left
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

  %% Saving                                                                      
  %                                                                                                                                             FUTURE PLOT SAVING
  processedImgfullpath = strcat(processedImgname,struct.MeasInfo{idxgral});
  % Explanation: saveas(variable,directory+filename,extension)
  saveas(gcf,strcat(processedImgfullpath),imgformat); % Saves the last shown figure
  
  % OLD:
  %   % Explanation: imwrite(variables,directory+filename+extension)
  %   imwrite(struct.expImgs{idxgral}, strcat(processedImgfullpath,dataformat));
  
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

%% Ask to leave figures open or not
answer = questdlg('Do you want to close all the figures?','Processing finished','yes','no');
if strcmp(answer,'yes') % Compare string
     close all;         
end
end
