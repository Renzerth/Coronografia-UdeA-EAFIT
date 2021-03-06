function [coorShiftX, coorShiftY] = f_centerMaskToSpot(referenceData, ...
         projectionMasks,mainDataCenter,mainDataRadius,monitorSize,TC,vid)
% A centering application designed to adjust optical vortex positioning.
% Based on a gradient step iterator of radial intensity the algorithm loops
% through the whole Screen in a coarse displacement looking for small
% intensity variations. Minor detections are used to stablish a scan
% area where displacements steps are refined and intensity variations are 
% analyzed in a finer dynamic manner until intensity behaves in a predicted
% vortex-like distribution around the spot axis.
%
%
% Inputs: referenceData (MxN) [uint8] Gray image of spot without vortex
%
%         projectionMask (PxL) [double] Digital hologram of a vortex
%         distribution to be presented at an SLM type screen
%
%         mainDataCenter (2x1) [double] Center coordinates of spot defined
%         in a previous procedure. See: findCircleShapedIntensity.m
%
%         mainDataRadius (0x1) [double] An estimative of the overall spot
%
%         monitorSize (2x1) [double]. Width and Height information of SLM
%         holds values for PxL size.
%
%         TC (0x1) [double] Topological Charge of the vortex mask
%
%         vid [videoinput] Video Input source object related to a feedback
%         camera.
%
% Outputs: shiftY (0x1) [double] Row Position of the vortex center related
%          to Y coordinates. Used for shifting the mask coordinates
%
%          shiftX (0x1) [double] Column Position of the vortex center
%          related to X coordinates. Used for shifting the mask coordinates
%
%
% Version 1.78 - for Matlab 2014b and onward distributions
%
% Author: Juan Jos� Cadavid Mu�oz
%         Engineering Physicist
%
% Contact: jcadav22@eafit.edu.co
%
% Affiliation: Universidad EAFIT 
%              Grupo �ptica Aplicada
%              Colombia - 2019.
%
% License: CC(4.0) BY-NC-SA
%
% Master in Applied Physics -- Optics Research
%% Input Check
adjustmentEnabled = false;
if length(projectionMasks) == 2
  projectionMask = projectionMasks{1}; % TC mask is used for selfcentering algorithm
  adjustmentEnabled = true;
end
%% Algorithm Settings
varBaseLine = 0.3;
criteriaTol = 0.5; % Ref: 0.5
histWeightFactor = 0.10;
spotFraction = 0.2; % Ref: 0.15 check
noiseValExcl = 15;
%% Set Reference
dataSize = size(referenceData);
[referenceRadialProfile] = f_getAverageRadialProfile(referenceData, ...
                                                  dataSize,mainDataCenter);
dataRange = noiseValExcl:TC*100/2;

%% Set Iterator
scanArea = round(spotFraction*mainDataRadius);
shiftStepX = f_getClosestMultiple(monitorSize(1),scanArea);
shiftStepY = f_getClosestMultiple(monitorSize(2),scanArea);
maxIterX = ceil(monitorSize(1)/shiftStepX);
maxIterY = ceil(monitorSize(2)/shiftStepY);
circIterator = @(maxIter,iter,n) (maxIter+1)*((iter+n)<0) + iter+n;

%% Figure-Screen properties
%%% Mask figure
offsetPixel = [1,1];
figureHandler = figure('Visible','off','MenuBar','none','Toolbar', ...
                       'none','Units','Pixels','color','black');
axesHandler = gca; % info:set handler property to 'figure' to revert 'none'
set(axesHandler,'Units','Pixels','Position',[offsetPixel monitorSize(1) ...
                                             monitorSize(2)]);
updateDisplayHandler = imagesc(zeros(fliplr(monitorSize))); axis off;
colormap('gray'); 
figureHandler.Visible = 'on'; axis on;
[~] = f_changeProjectionMonitor('Restore');

%%% Reference radial profile figure
FigProfHandler = figure('color','white'); 
analysisPlotHandler = plot(referenceRadialProfile);ylim auto; grid on;
% This title applies since later on the relative difference will be
% dynamically plotted here

axis square; title(['Difference between reference (tc=0) and tc=',num2str(TC)]);
xlabel('Radial distance from the center of the spot [pixels]');
ylabel('Relative difference of intensities [a.u.]')

%% Specific region of scanning
% rectInfo = imrect(axesHandler);
% position = round(rectInfo.getPosition);
% shiftStepX = getClosestMultiple(position(3), ...
%                                 round(spotFraction*position(3)));
% shiftStepY = getClosestMultiple(position(4), ...
%                                 round(spotFraction*position(4)));
% maxIterX = position(3)/shiftStepX;
% maxIterY = position(4)/shiftStepY;
% searchOffsetX = position(1);
% searchOffsetY = position(2);

%% Centering logic properties
searchOffsetX = 0;
searchOffsetY = 0;
backwardIter = 2;
dynamicEval = 0;
scanAreaFactor = 2*backwardIter;
previousVal = 0;
refiningIterations = 3; % To large will lead to infinite state -> 4
updateEnabled = false;

%% Centering logic loop
for times = 1: refiningIterations
  for iterationsY = 0:maxIterY
    for iterationsX = 0:maxIterX
      pause(0.05);
      currentFrame = getsnapshot(vid);
      wait(vid);
      [criteriaValue, relativeChange, ~] = f_getDistMetrics(currentFrame,dataSize, mainDataCenter, referenceRadialProfile, dataRange);
      valComparison = (criteriaValue - previousVal);
      previousVal = criteriaValue;
      disp([criteriaValue dynamicEval]); % criteriaValue: average in the vortex region of the relative difference of the profile
                                                           % dynamicEval: defines when there's anintensity change
      
      if abs(valComparison) >= 0.1 || criteriaValue >= criteriaTol
        disp('Spot rapid change - Check')
        break;
      end
      
      localY = iterationsY*shiftStepY + searchOffsetY; % If breaks before updating, it holds the variation coordinates (n-1)
      localX = iterationsX*shiftStepX + searchOffsetX;
      set(updateDisplayHandler,'CData',circshift(projectionMask,[localY, localX]));
      set(analysisPlotHandler,'YData',relativeChange);
    end
        dynamicEval = (varBaseLine + abs(times/10 - histWeightFactor*(criteriaValue + valComparison)) + 0.1*(1-exp(-(2^(times-1)/TC))));
    if criteriaValue >= dynamicEval
      disp('Spot high change - stop')
      updateEnabled = true;
      break;
    end
  end
  
  if updateEnabled
    searchOffsetY = circIterator(maxIterY,iterationsY,-(backwardIter - 1))*shiftStepY + searchOffsetY;
    searchOffsetX = circIterator(maxIterX,iterationsX,-backwardIter)*shiftStepX + searchOffsetX;
%     searchOffsetY = searchOffsetY + monitorSize(1)*(searchOffsetY == 0 && times > 1); % Hysteresis Loss Correction
%     searchOffsetX = searchOffsetX + monitorSize(2)*(searchOffsetX == 0 && times > 1);
    maxIterX = scanAreaFactor*shiftStepX;
    maxIterY = scanAreaFactor*shiftStepY;
    shiftStepX = round(shiftStepX/(2^(times+1)));
    shiftStepY = round(shiftStepY/(2^(times+1)));
    maxIterX = ceil(maxIterX/shiftStepX);
    maxIterY = ceil(maxIterY/shiftStepY);
    iterationsX = 0;
    iterationsY = 0;
    updateEnabled = false;
    pause(1); % Wait for intensity to be properly registered after scan area change
  else
    if ~adjustmentEnabled
      error('Failed to find representative radial changes during iteration.');
    end
  end
end
%% Manual adjustment of the mask
manualShiftXcoord = 0;
manualShiftYcoord = 0;
if adjustmentEnabled
  % Instruction for the user
  disp(['Move the sliders to adjust the mask"s position and then close'...
     ' the mask window when done']);
  set(updateDisplayHandler,'CData',circshift(projectionMasks{2},[localY, localX]));
  [manualShiftYcoord, manualShiftXcoord] = f_adjustSLMPositioning( ...
  figureHandler,updateDisplayHandler,analysisPlotHandler,vid,dataSize, ...
  mainDataCenter,referenceRadialProfile,dataRange);
end

%% Coordinates correction
circShiftX = localX + manualShiftXcoord;
circShiftY = localY + manualShiftYcoord;
[coorShiftX, coorShiftY] = f_circShiftToCart(circShiftX,circShiftY,monitorSize(1),monitorSize(2)); % Referred % Overall Coordinates from screen origin

%% Close the profile figure
if ishandle(FigProfHandler)
  close(FigProfHandler);
end

end