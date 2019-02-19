%% Projection Gray Level ranges
minGrayDepth = 0;
maxGrayDepth = 10;
levelShift = 0;
%% Screen Coordinates
screenIndex = 3;
[X,Y,~,monitorSize] = makeScreenCoordinates(screenIndex);
%% Get Center position
% [~,relativeCoord] = TOOLS.computeCenterFromImageMin();
% shiftCorr = [0.113,0.018];
%% Center Shift
% shiftX = 0.5 - relativeCoord(2) + shiftCorr(2);
% shiftY = 0.5 - relativeCoord(1) + shiftCorr(1);
shiftX = 0;
shiftY = 0;
%% Vortex Phase
TC = 10;
maxGrayLevel = 255;
% [X,Y] = meshgrid(-1:2/255:1); % Is not scaled
[spiralPhase,~] = cart2pol(X-shiftX,Y-shiftY);
projectionMask = angle(exp(1i*spiralPhase*TC));
grayLevels = -pi:2*pi/(maxGrayLevel):pi;
projectionMask = discretizeMask(projectionMask,grayLevels);
%% To uint 8 Values Scaling
projectionMask = scaleMatrix(projectionMask,minGrayDepth,maxGrayDepth) + levelShift;
%% Figure-Screen properties
close(gcf);
offsetPixel = [1,1];
figureHandlerA = figure('Visible','off','MenuBar','none','Toolbar','none');
figureHandlerA.Units = 'Pixels';
set(gca,'Units','Pixels');
set(gca,'Position',[offsetPixel monitorSize(1) monitorSize(2)]);
%% Project Object
image(projectionMask); axis off; colormap(gray(maxGrayDepth));
figureHandlerA.Visible = 'on';
%% Restore default figure
[~] = changeProjectionMonitor('Restore');