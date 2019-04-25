function [x,y,Hprof,Vprof,maxX,maxY] = f_makeImageProfile(x,y,midX,midY, ...
    dataArray,shiftCart, oneSideProfile, varargin)
% returns the profiles of a 2D data array (or matrix), i.e., an image or a
% 2D map. These profiles are horizontal and vertical
%
% Inputs:
%  x,y: cartesian coordinates for dataArray
% midX:
% midY:
%  dataArray: 2D numerical array (or matrix): image/mask for example
%  size of the profile
%  shiftCart: should only be in [-0.99 , 0.99]
%
% Optional Inputs (Enables plotting)
%  tit: title for dataArray
%  xlab,ylab: x and y labels (strings)
%  plotData: plots dataArray and lines of the profiles
%  plotH,plotV: booleans to plot the horizontal and the vertical profiles
%  oneSideProfile: 0 (no); 1(1st and 4th quadrants); 2(2nd and 3rd
%  quadrants). For cases 1 and 2, the origin is always on the center of the
%  image (determined by shiftCart)
%  tol: discrete index tolerance in order to symmetrically dicrease the
%
% Outputs:
%  [Hprof,Vprof]: array of the profiles
%% Program settings and input check
if nargin == 14 && prod(cellfun(@(inputs) ischar(inputs), varargin{1:3})) && prod(isnumeric(varargin(4:7)))
    tit = varargin{1};
    xlab = varargin{2};
    ylab = varargin{3};
    plotData = varargin{4};
    plotH = varargin{5};
    plotV = varargin{6};
    tol = varargin{7};
    plottingEnabled = true;
else
    plottingEnabled = false;
end

%% Max points of the array (with the shift application)
% Takes into account an even/odd size of the size of dataArray
[maxY, maxX] = size(dataArray); % Size of the array

%% Mid points of the array
% midX = round((maxX+1)/2) + dcShift*mod(maxX,2); % x mid point
% midY = round((maxY+1)/2) + dcShift*mod(maxY,2); % y mid point

%% Shift of the mask scaling
shiftX = shiftCart(2); % Cartesian shift in x
shiftY = shiftCart(1); % Cartesian shift in y
% shiftX = round(midX*shiftX); % Percentage of the half support of the mask
% shiftY = round(midY*shiftY); % Percentage of the half support of the mask

%% Mid points of the array (with the shift application)
% The extreme parts of the profile are conserved and the midpoints shift
% The signs here are the opposite of the Shift application section of the
% function "f_DefineSpace" since they compensate this shift
midX = midX + shiftX; % Shifted X for the SLM
midY = midY - shiftY; % Shifted Y for the SLM

%% Profile coordinates
Hx = [1,maxX]; % Horizontal x components: (xi,xf)
Vy = [1,maxY]; % Vertical y components: (yi,yf)
Hy = [midY,midY]; % Horizontal y components: (yi,yf)
Vx = [midX,midX]; % Vertical x components: (xi,xf)

%% Calculation of the profiles
Hprof = improfile(dataArray,Hx,Hy); % Horizontal profile
Vprof = improfile(dataArray,Vx,Vy); % Vertical profile

%% One-side profile extraction
switch oneSideProfile
    case 0
        % No change: full-sided profile
    case 1 %  Profiles in the 1st and 4th quadrants
        Hprof = Hprof(midX:end); % From the middle to the end
        Vprof = Vprof(midY:end); % From the middle to the end
        x = x(midX:end); % From the middle to the end
        y = y(midY:end); % From the middle to the end
    case 2 % Profiles in the 2nd and 3rd quadrants
        Hprof = fliplr(Hprof(1:midX)'); % From the middle to the end
        Vprof = fliplr(Vprof(1:midY)'); % From the middle to the end
        x = x(1:midX); % From the middle to the end
        y = y(1:midY); % From the middle to the end
end

%% Plotting
if plottingEnabled
    f_plotLinearProfiles(dataArray,x,y,tit,xlab,ylab,plotData,plotH,plotV,tol)
end
end