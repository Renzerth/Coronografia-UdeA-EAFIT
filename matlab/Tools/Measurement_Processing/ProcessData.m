% [A] = f_DataProcessing(imgfullpath,savetype,pathSep,dataformat,cameraPlane,n,PP,M,f,ProcessedDir)
%% Post-processing of the data (application of the metric of the degree of
%%% extintion)

%% Measurements initializing
% Copyright PhD student Jens de Pelsmaeker VUB B-PHOT 2018,Brussels,Belgium
pause(wait) % Seconds before measuring as a safety measurement
t1_dt = datetime; % store time
disp('Processing started:'); disp(t1_dt)
processedImgname = strcat(ProcessedDir,pathSep,'processed_', ...
  cameraPlane,'_');

% for idxgral = 1:totalImgs
%% Loading
% load(directory+filename,variables)

A = load(imgfullpath,'expImgs'); % imgfullpath comes from AutomatMeasure

%% Cartesian coordinates with the lambda/D scaling
[ySize, xSize] = size(distribution);
xpix = 1:xSize;
ypix = 1:ySize;
xang = f_scalePix2DiffAng(xpix,n,PP,M,f);
yang = f_scalePix2DiffAng(ypix,n,PP,M,f);

%% Processsing of the image
switch metricSel
  case 1
    tit = 'Encircled Energy Factor metric';
    [energy,radialIntensity] = f_calculateEEF(angle(mask),n,PP,M,f,shiftCart,metricProfile,tit);
end

%% Saving
processedImgfullpath = strcat(processedImgname,MeasInfo{idxgral});
%  if savetype == 1 % .mat format
%   % save(directory+filename,variables)
%   save(processedImgfullpath,'A'); % .mat
%  else % savetype = 2. dataformat is used
%   % imwrite(variables,directory+filename+extension)
% imwrite(expImgs{idxgral}, strcat(processedImgfullpath,dataformat));
%  end
%end

%% End of the processing
% Author: PhD student Jens de Pelsmaeker VUB B-PHOT 2018, Brussels, Belgium
% MATLAB built in:
t2_dt = datetime;
disp('Processing finished:'); disp(t2_dt)
time = t2_dt - t1_dt;
disp('Processing took: '); % datestr(time,'SS') ' seconds'])
disp(time)
% end