function [X,scaledY,R,monitorSize] = f_makeScreenCoordinates(screenIndex,varargin)
% Inputs:
% screenIndex: screen number selector
% varargin: function activated (1) or not (0)

%% Input Verification
if nargin == 2 % number of inputs
    if isa(varargin{1},'logical')
        relativeCoordSelect = varargin{1};
    else
        relativeCoordSelect = true;
    end
else
    relativeCoordSelect = true;
end

%% Projection monitor properties
[monitorSize,~] = f_changeProjectionMonitor(screenIndex);
horizontalHalfSize = ceil(monitorSize(1)/2);
verticalHalfSize = ceil(monitorSize(2)/2);
aspectRatio = monitorSize(1)/monitorSize(2);

%% Space definition
if relativeCoordSelect
    horSpaceArray = -1:2/(monitorSize(1)-1):1;
    verSpaceArray = -1:2/(monitorSize(2)-1):1;
else
    horSpaceArray = -horizontalHalfSize:horizontalHalfSize-1;
    verSpaceArray = -verticalHalfSize:verticalHalfSize-1;
end
[X,Y] = meshgrid(horSpaceArray,verSpaceArray);
scaledY = Y/aspectRatio;
R = sqrt(X.^2 + scaledY.^2);
end