function [X,Y,AspectRatio,monitorSize] = f_MakeScreenCoords(scrnIdx,enablechange,varargin)
% Inputs:
% screenIndex: screen number selector. In [1,N] with N the # of screens
% varargin: relativeCoordSelect. false: pixels space; true: unitary spaces
%
% Outputs:
% X and Y vectors of the size of the screem
% AspectRatio: of the screen
% monitorSize: of the screen

%% Input Verification
if nargin == 3 % number of inputs
    if isa(varargin{1},'logical')
        relativeCoordSelect = varargin{1};
    else
        relativeCoordSelect = true; % Unitary space
    end
else
    relativeCoordSelect = true; % Unitary space
end

%% Projection monitor properties
[monitorSize,~] = f_changeProjectionMonitor(scrnIdx,enablechange);
horizontalHalfSize = ceil(monitorSize(1)/2);
verticalHalfSize = ceil(monitorSize(2)/2);
AspectRatio = monitorSize(1)/monitorSize(2); % Horizontal over vertical
% The aspect ratio of a rectangle is the ratio of its longer side to its
% shorter side – the ratio of width to height,[1] when the rectangle is 
% oriented as a "landscape". 

%% Space definition
if relativeCoordSelect % Unitary space
    horSpaceArray = -1:2/(monitorSize(1)-1):1;
    verSpaceArray = -1:2/(monitorSize(2)-1):1;
else % Pixels
    horSpaceArray = -horizontalHalfSize:horizontalHalfSize-1;
    verSpaceArray = -verticalHalfSize:verticalHalfSize-1;
end
[X,Y] = meshgrid(horSpaceArray,verSpaceArray);
end