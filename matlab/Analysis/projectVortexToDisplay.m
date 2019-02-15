% Questions:
% Differente between shiftX,Y and offsetPixel
% levelShift seems not to be linear

%% Projection Gray Level ranges
minGrayDepth = 0; % Ref: 0
maxGrayDepth = 100; % Ref: 255
levelShift = 0; % Ref: 0. Seems to be non-linear or better not to use it

%% Screen Coordinates
screenIndex = 1; % Number of screens connected to the pc
[X,Y,~,monitorSize] = f_makeScreenCoordinates(screenIndex);

%% Get Center position
% [~,relativeCoord] = f_computeCenterFromImageMin();
% shiftCorr = [0.113,0.018];

%% Center Shift
% shiftX = 0.5 - relativeCoord(2) + shiftCorr(2);
% shiftY = 0.5 - relativeCoord(1) + shiftCorr(1);
shiftX = 0;
shiftY = 0;

%% Vortex Phase
k = 10; % Bits for grey levels
k = 2^k; % Resolution of the grey levels
tc = 1; % Topological charge
maxGrayLevel = 100; % Adjusts the # of gray levels of the mask. Must be
                    % smaller than (maxGrayDepth-minGrayDepth)
[X,Y] = meshgrid(-1:2/(k-1):1);
[spiralPhase,~] = cart2pol(X-shiftX,Y-shiftY);
projectionMask = angle(exp(1i*spiralPhase*tc));
grayLevels = -pi:2*pi/(maxGrayLevel):pi;
projectionMask = f_discretizeMask(projectionMask,grayLevels);

%% To uint 8 Values Scaling
projectionMask = f_scaleMatrix(projectionMask,minGrayDepth,maxGrayDepth) + levelShift;

%% Figure-Screen properties
close(gcf);
offsetPixel = [0,0];
figureHandlerA = figure('Visible','off','MenuBar','none','Toolbar','none');
figureHandlerA.Units = 'Pixels';
set(gca,'Units','Pixels');
set(gca,'Position',[offsetPixel monitorSize(1) monitorSize(2)]);

%% Project Object
image(projectionMask); axis off; colormap(gray(maxGrayDepth));
figureHandlerA.Visible = 'on';

%% Restore default figure
[~] = f_changeProjectionMonitor('Restore');