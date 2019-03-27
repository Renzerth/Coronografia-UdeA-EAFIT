function f_ProcessData(measfullpath,refmeasfullpath,ProcessedDir,pathSep,infoDelim, ...
dataformat,cameraPlane,totalImgs,AiryFactor,metricSel,metricProfile, ...
shiftCart,beepSound,L,NA)
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
I = load(measfullpath); % Loads all the measured images and their info
                        % Two variables are loaded in the structure: 
                        % "expImgs" and "MeasInfo"
                        
%%  Loading the reference measurement
refmeas = imread(refmeasfullpath);          
refmeas = im2double(refmeas);

%% Find center of the PSF or the Lyot image

                        
%% Cartesian coordinates with pixel units
[ySize, xSize] = size(I.expImgs{1}); % All images assumed of the same size
% xpix = 1:xSize; % Pixels start in 1
% ypix = 1:ySize; % Pixels start in 1

%% Airy radius calculation
airyBool = 1;
switch airyBool
    case 1 % Measured from the reference image
       AiryDiskPixX = 123; % Just an example
       AiryDiskPixY= 79; % Just an example                                  % ORG WITH TEMPORAL.m
       % Leer, binarizar, erode
       
%        AiryDiskSpatialX = AiryDiskPixX*PP; % PP: camera's pixel pitch
%        AiryDiskSpatialY = AiryDiskPixY*PP;
    case 2 % Theoretical: diffraction-limited systems
       % NA: numerical aperture of the lens
        % L: wavelength in um
        AiryDiskSpatial = 0.61*L/NA; % um
        AiryDiskSpatialX = AiryDiskSpatial;
        AiryDiskSpatialY = AiryDiskSpatial;
        
        AiryDiskPixX = AiryDiskSpatialX/PP; % PP: camera's pixel pitch
        AiryDiskPixY = AiryDiskSpatialY/PP;
        
        % In reality, there are aberrations and the real radius can be of
        % about 60 times the theoretical one
    case 3 % From the EEC factor
        % When it is the 70%
end

%% Pixel airy radius
AiryMultiplicityX = xSize/AiryDiskPixX; % NUmber of airy disks 
AiryMultiplicityY = ySize/AiryDiskPixY;

xpix = -xSize/2 : 1 : xSize/2 - 1; % pixels
ypix= -ySize/2 : 1 : ySize/2 - 1; % pixels                   % MAYBE USE MIDX,MIDY

xangL_D = xpix/(AiryMultiplicityX*8);
yangL_D =  ypix/(AiryMultiplicityY*8);

%% Cartesian coordinates with the lambda/D scaling (diffraction angle)
% xangL_D = f_scalePix2DiffAng(xpix,AiryFactor);  
% yangL_D = f_scalePix2DiffAng(ypix,AiryFactor);

%% Cartesian coordinates with the arcsecond scaling (diffraction angle)
xangArcs = f_LambdaDToarcsec(xangL_D);
yangArcs = f_LambdaDToarcsec(yangL_D);

%% Metric-specific default parameters
xlab = 'Angular position [\lambda/D]';
ylab = 'Angular position [\lambda/D]';
plotData = 1; % Shows the profile lines. Ref: 1
plotH = 0; % Not needed for the metric. Ref: 0
plotV = 0; % Not needed for the metric. Ref: 0
oneSideProfile = 2; % Specifically needed for this metric. Ref: 1
dcShift = 0; % Only used for spectra (Fourier analysis)
tol = 0; % 0: no need to symmetrically truncate the profile. Ref: 0

for idxgral = 1:totalImgs
  %% Processsing of the image
  
  % ACTUALLY, CALCULATE ALL THE METRICS BUT ONLY PLOT THE WANTED ONES
  
  switch metricSel
    case 1 % Throughput: Encircled Energy Factor metric
      % Camera: PSF.
      tit = 'Encircled Energy Distribution of Intensity';
      [~,~] = f_calculateEEF(xangL_D,yangL_D,I.expImgs{idxgral}, ...
      shiftCart,metricProfile,tit,xlab,ylab,plotData,plotH,plotV, ...
      oneSideProfile,dcShift,tol);
    case 2 % Throughput gradient
      tit = 'Throughput gradient'; 
      % gradient [returns n elements] or diff [returns n-1 elements]
       
    case 3 % Power suppresion in the airy disk
      tit ='Power suppresion in the airy disk';    
      % Needs case 1 as input
      
    case 4 % SNR
      tit = 'Signal-to-Noise Ratio';
      [~] = f_calculateSNR(xangL_D,yangL_D,I.expImgs{idxgral},refmeas,shiftCart, ...
      metricProfile,tit,xlab,ylab,plotData,plotH,plotV,oneSideProfile, ...
      dcShift,tol);
    case 5 % MSE
      tit = 'Mean Squared Error';
        
    otherwise
      error('Select a valid metric');
  end

  %% Saving
  processedImgfullpath = strcat(processedImgname,I.MeasInfo{idxgral});
  % Explanation: imwrite(variables,directory+filename+extension)
  imwrite(I.expImgs{idxgral}, strcat(processedImgfullpath,dataformat));
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
