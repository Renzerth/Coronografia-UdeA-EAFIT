function f_ProcessData(measfullpath,refmeasfullpath,ProcessedDir,pathSep,infoDelim, ...
dataformat,cameraPlane,totalImgs,AiryFactor,metricSel,metricProfile, ...
shiftCart,beepSound,L,NA,PP)
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
[ySize, xSize] = size(refmeas); % All images assumed of the same size as the refmeas
% xpix = 1:xSize; % Pixels start in 1
% ypix = 1:ySize; % Pixels start in 1

%% Theoretical: diffraction-limited systems
% NA: numerical aperture of the lens
% L: wavelength in um
AiryDiskSpatial = 0.61*L/NA; % um
AiryDiskSpatialX = AiryDiskSpatial;
AiryDiskSpatialY = AiryDiskSpatial;
disp(strcat("Theoretical Airy's radius in um: ", num2str(AiryDiskSpatialX), ' (x) and ', num2str(AiryDiskSpatialY), ' (y)'));

AiryDiskPixX = AiryDiskSpatialX/PP; % PP: camera's pixel pitch
AiryDiskPixY = AiryDiskSpatialY/PP;
disp(strcat("Theoretical Airy's radius in pix: ", num2str(AiryDiskPixX), " (x) and ", num2str(AiryDiskPixY), ' (y)'));

% In reality, there are aberrations and the real radius can be of
% about 60 times the theoretical one

%% Airy radius measured from the reference image
AiryDiskPixX = aproxRadius; % Just an example
AiryDiskPixY= aproxRadius; % Just an example
AiryDiskSpatialX = AiryDiskPixX*PP; % PP: camera's pixel pitch
AiryDiskSpatialY = AiryDiskPixY*PP;
disp(strcat("Estimated Airy's radius in um: ", num2str(AiryDiskSpatialX), " (x) and ", num2str(AiryDiskSpatialY), ' (y)'));
disp(strcat("Estimated Airy's radius in pix: ", num2str(AiryDiskPixX), " (x) and ", num2str(AiryDiskPixY), ' (y)'));

%% Airy radius from the EEC factor
% When it is the 70%

%% Pixel airy radius
% AiryMultiplicityX = xSize/AiryDiskPixX; % Number of airy disks 
% AiryMultiplicityY = ySize/AiryDiskPixY;

%% Symmetric pixels
xpix = -xSize/2 : 1 : xSize/2 - 1; % pixels
ypix= -ySize/2 : 1 : ySize/2 - 1; % pixels                   % MAYBE USE MIDX,MIDY

% xangL_D = xpix/(AiryMultiplicityX*8);
% yangL_D =  ypix/(AiryMultiplicityY*8);

%% Cartesian coordinates with the lambda/D scaling (diffraction angle)
% xangL_D = f_scalePix2DiffAng(xpix,AiryFactor);  
% yangL_D = f_scalePix2DiffAng(ypix,AiryFactor);

%% Cartesian coordinates with the arcsecond scaling (diffraction angle)
% xangArcs = f_LambdaDToarcsec(xangL_D);
% yangArcs = f_LambdaDToarcsec(yangL_D);

%% Lyot's spot size (main radius)
apRad = 2.54; % Aperture radius [cm]

Lyotimg = imread('/home/labfisica/Dropbox/Coronógrafo_2018-1_Samuel/6_Photos/Project development/5-Measurement-Results/18_data_ref-self_centering/data_ref_2.bmp');
Lyotimg = rgb2gray(Lyotimg);

PSFimg = imread('/home/labfisica/Dropbox/Coronógrafo_2018-1_Samuel/two focal planes.bmp');


[ySize, xSize] = size(PSFimg); 
% xpix = -xSize/2 : 1 : xSize/2 - 1; % pixels
xpix = 1:xSize;
ypix= 1:ySize;
% ypix= -ySize/2 : 1 : ySize/2 - 1; % pixels   

drawing = false;
PP = PP*1e-6; % um to m
apRad = apRad*1e-2; % cm to m
[~,mainLyotRadius,~] = findCircleShapedIntensity(Lyotimg,drawing);
[apRadpix] = computePupilPixelSize(mainLyotRadius, PP,apRad);

%%  lambda/D factor falco-matlab reference
% It is scalled with respect to the jinc zeros
[ySize,xSize] = size(Lyotimg);
NpadX = xSize; % Camera's x pixel size
NpadY = ySize; % Camera's y pixel size
xlamOverD = NpadX/(2*apRadpix);
ylamOverD = NpadY/(2*apRadpix);
xvals_fp = -NpadX/2:NpadX/2-1;
x_L_D = xvals_fp/xlamOverD;
yvals_fp = -NpadY/2:NpadY/2-1;
y_L_D = yvals_fp/ylamOverD;

%% Find center of the PSF image (binarization)
[xangL_D,yangL_D,regionCentroid,aproxRadius] = f_approximateSpotSize(PSFimg);

%% Profile
plotData = 1; % Shows the profile lines. Ref: 1
plotH = 0; % Not needed for the metric. Ref: 0
plotV = 0; % Not needed for the metric. Ref: 0
oneSideProfile = 1; % Specifically needed for this metric. Ref: 1
dcShift = 0; % Only used for spectra (Fourier analysis)
tol = 0; % 0: no need to symmetrically truncate the profile. Ref: 0
tit = 'tit';
xlab = 'Angular position [\lambda/D]';
ylab = 'Angular position [\lambda/D]';

midX = round((xSize+1)/2); % x mid point
midY = round((ySize+1)/2); % y mid point

regionCentroid = round(regionCentroid);

sX = regionCentroid(2);
sY = regionCentroid(1);
 
sX =sX- midX;
sY = midY-sY;
 
regionCentroid = [sY sX];

regionCentroid =regionCentroid./(2*[ySize,xSize]); % Nearest pixel

[x,y,Hprof,Vprof,~,~,~,~] = f_makeImageProfile(xpix,ypix,PSFimg,...
tol,regionCentroid,tit,xlab,ylab,plotData,plotH,plotV,oneSideProfile,dcShift);




%% Metric-specific default parameters
xlab = 'Angular position [\lambda/D]';
ylab = 'Angular position [\lambda/D]';
plotData = 1; % Shows the profile lines. Ref: 1
plotH = 0; % Not needed for the metric. Ref: 0
plotV = 0; % Not needed for the metric. Ref: 0
oneSideProfile = 1; % Specifically needed for this metric. Ref: 1
dcShift = 0; % Only used for spectra (Fourier analysis)
tol = 0; % 0: no need to symmetrically truncate the profile. Ref: 0

for idxgral = 1:totalImgs
  %% Processsing of the image
  
  % ACTUALLY, CALCULATE ALL THE METRICS BUT ONLY PLOT THE WANTED ONES
  
  switch metricSel
    case 1 
     %% Throughput: Encircled Energy Factor metric
      % Camera: PSF.
      tit = 'Encircled Energy Distribution of Intensity';
      [energy,radialIntensity,cartcoord,titprof] = f_calculateEEF(xpix,ypix,struct.expImgs{idxgral}, ...
      shiftCart,metricProfile,tit,xlab,ylab,plotData,plotH,plotV, ...
      oneSideProfile,dcShift,tol);
     
     %% Throughput gradient
      tit = 'Throughput gradient'; 
      % gradient [returns n elements] or diff [returns n-1 elements]
      normIntensity = radialIntensity./max(radialIntensity);
      GradEnergy = gradient(normIntensity);
      
     %% Plot of the gradient of the EEF and its corresponding intensity pattern
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
      % Needs case 1 as input
      
    case 2 % SNR
      tit = 'Signal-to-Noise Ratio';
      [~] = f_calculateSNR(xpix,ypix,struct.expImgs{idxgral},refmeas,shiftCart, ...
      metricProfile,tit,xlab,ylab,plotData,plotH,plotV,oneSideProfile, ...
      dcShift,tol);
    case 3 % MSE
      tit = 'Mean Squared Error';
        
    otherwise
      error('Select a valid metric');
  end

  %% Saving
  processedImgfullpath = strcat(processedImgname,struct.MeasInfo{idxgral});
  % Explanation: imwrite(variables,directory+filename+extension)
  imwrite(struct.expImgs{idxgral}, strcat(processedImgfullpath,dataformat));
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
N = 4; % Number of beeps
f_EndBeeps(N,beepSound);

%% Ask to leave figures open or not
answer = questdlg('Do you want to close all the figures?','Processing finished','yes','no','no');
if strcmp(answer,'yes') % Compare string
     close all;         
end
end
