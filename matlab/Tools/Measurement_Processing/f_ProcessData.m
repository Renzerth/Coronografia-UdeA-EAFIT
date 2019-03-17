function f_ProcessData(measfullpath,ProcessedDir,pathSep,infoDelim, ...
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

%% Loading
% Explanation: load(directory+filename,variables)
I = load(measfullpath); % Loads all the measured images and their info
                        % Two variables are loaded in the structure: 
                        % "expImgs" and "MeasInfo"

%% Cartesian coordinates with the lambda/D scaling (diffraction angle)
[ySize, xSize] = size(I.expImgs{1}); % All images assumed of the same size
xpix = 1:xSize;
ypix = 1:ySize;
xangL_D = f_scalePix2DiffAng(xpix,AiryFactor);
yang_D = f_scalePix2DiffAng(ypix,AiryFactor);

%% Cartesian coordinates with the arcsecond scaling (diffraction angle)
xangArcs = f_LambdaDToarcsec(xangL_D);
yangArcs = f_LambdaDToarcsec(yang_D);

for idxgral = 1:totalImgs
  %% Processsing of the image
  switch metricSel
    case 1
      tit = 'Encircled Energy Factor metric';
      [energy,radialIntensity] = f_calculateEEF(I.MeasInfo{idxgral},shiftCart, ...
                                                metricProfile,tit);
  end

  %% Saving
  processedImgfullpath = strcat(processedImgname,I.expImgs{idxgral});
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
disp(time)
end
