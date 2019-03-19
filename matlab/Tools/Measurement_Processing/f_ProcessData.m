function f_ProcessData(measfullpath,refmeasfullpath,ProcessedDir,pathSep,infoDelim, ...
dataformat,cameraPlane,totalImgs,AiryFactor,metricSel,metricProfile, ...
shiftCart)
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
                        
%% Cartesian coordinates with pixel units
[ySize, xSize] = size(I.expImgs{1}); % All images assumed of the same size
xpix = 1:xSize; % Pixels start in 1
ypix = 1:ySize; % Pixels start in 1

%% Airy radius
airyBool = 1;
switch airyBool
    case 1 % Theoretical
       Rairy = 1;
    case 2 % Measured from the reference image
        NA = 0.1;
        L = 0.6328;
        Rairy = 1.22*L/NA;
end

%% Cartesian coordinates with the lambda/D scaling (diffraction angle)
xangL_D = f_scalePix2DiffAng(xpix,AiryFactor);  
yangL_D = f_scalePix2DiffAng(ypix,AiryFactor);

%% Cartesian coordinates with the arcsecond scaling (diffraction angle)
xangArcs = f_LambdaDToarcsec(xangL_D);
yangArcs = f_LambdaDToarcsec(yangL_D);

for idxgral = 1:totalImgs
  %% Processsing of the image
  switch metricSel
    case 1 % Throughput: Encircled Energy Factor metric
      % Camera: PSF.
      tit = 'Encircled Energy Distribution of Intensity';
      [~,~] = f_calculateEEF(xangL_D,yangL_D,I.expImgs{idxgral}, ...
                             shiftCart,metricProfile,tit);
      case 2 % Throughput gradient
        tit = 'Throughput gradient';  
      case 3 % Power suppresion in the airy disk
        tit ='Power suppresion in the airy disk';    
      case 4 % SNR
        tit = 'Signal-to-Noise Ratio';
        % refmeas
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
end
