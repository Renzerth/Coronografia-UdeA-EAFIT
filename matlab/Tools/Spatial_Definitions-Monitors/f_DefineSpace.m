function [rSize,x,y,Xslm,Yslm,rSLM,phiSLM,Xpc,Ypc,rPC,phiPC,ZernikeSize, ...
monitorSize,coordType] = f_DefineSpace(sSupport,sSize,shiftCart, ...
shiftBool,pixSize,scrnIdx,circularMask,coordType,MaxMask,plotMask,maskSel)
% Inputs:
%  sSupport: full side-length of the SLMs (or unitary without an SLM)
%  sSize: bits for grey levels; 2^k is the resolution (size of x and y)
%     Default: 10. Size is calculated as 2^k - 1
%     Only works when coordType = 1
%  shiftCart: [yshift,xshift], works when shiftBool = 1
%             Percentages of movement of the total size of the mask 
%             (cartesian coordinates convention). Calibrated with: s = +1;
%             ph0 = 0, tc = 1. Ranges per shift: [0,100] (percentage)  
%  shiftBool: only shifts when plotMask = 2
%             0: shift deactivated [for exporting masks]
%             1: shift activated [SLM displaying]
%             2: self-centering algorithm
%  pixSize: SLM pixel's size in um
%  scrnIdx: screen number selector. In [1,N] with N the # of screen
%  circularMask: Only works when coordType = 2
%           0: The mask presents an elliptical form when in the full screen
%           1: The mask presents a circular form when in the full screen
%           On both cases full screen means that plotMask = 2
%           It is always applied for Zernike masks (maskSel=5,6) either for 
%           the PC or for the SLM
%  coordType: type of calculation of the spatial coordinates. def: 2 
%    -1: size defined by the user, space support defined by the SLM to use
%    -2: size defined by the resolution of the selected screen
%  MaxMask: defines if the mask should be maximized when coordType = 1
%           -0: custom-size mask that depends on the variable sSize   
%           -1: maximizes the mask for coordType = 1
%           -2: maximized mask but keeping its rectangular fashion
%  plotMask:  no (0); on the screen (1); on the SLM (2); on the screen, but
%             a surface (3)
%  maskSel: selects a specific mask
%  coordType: input and output since it can change for maskSel = 5 and 6
%
% Outputs:
% rSize: radius for the circular (or elliptical) pupil truncation
% x,y: cartesian coordinates
% Xslm,Yslm,rSLM,phiSLM,Xpc,Ypc,rPC,phiPC: meshgrids of 2D cartesian and 
%                                          polar coordinates of both the 
%                                          SLM and the PC
% slm: shifts possibly applied; pc: shifts are never applied
% ZernikeSize: screen size for the Zernike polynomials generation
% monitorSize: resolution of the monitor selected by the scrnIdx

%% Screen coordinates and aspect ratio
  enablechange = false; 
  % false: won't change default figure display monitor. Leave this value as
  % zero always as the figure display monitor will be changed later on.
  % This is changed if one wants to display in the SLMs when plotMask = 2
  [Xcoord2,Ycoord2,AspectRatio,monitorSize] = ...
  f_MakeScreenCoords(scrnIdx,enablechange); % Calculates the monitor size

%% Zernike mask changes coordType from 2 to 1
% In order to apply circularMask for Zernike, coordType 1 is used with the
% smallest size of the monitor
if coordType == 2 && circularMask == 1 && (maskSel == 5 || maskSel == 6)
       sSize = min(monitorSize); % Minimum since that's the size of a 
                                 % square that can fit inside the screen's 
                                 % rectangle
       coordType = 1; % sSize will be used in this coordinate type
end

%% Coordinate type selection  
switch coordType
 case 1 % Size defined by the user, space support defined by the SLM to use
  % This applies when one won't project a full screen mask and a desired 
  % image resolution is wanted 
  %% Spatial definitions (size is user defined)
  sSupport = sSupport/2; % Half support of the SLM window in cm
  spaceVector = -sSupport:2*sSupport/(sSize-1):sSupport; 
  % Symmetric space
  [Xcoord1,Ycoord1] = meshgrid(spaceVector); % A symmetric grid: 2D 
                                             % Cartesian coordinates
  x = spaceVector; % Cartesian x-vector
  y = x; % Cartesian y-vector: square grid
  % Xslm = X;
  monitorSize = size(Xcoord1); % Square coordinates: [sSize,sSize]
  X = Xcoord1; Y = Ycoord1;
  
 case 2 % Size defined by the resolution of the selected screen     
  % This applies when a full screen mask will be displayed for the SLM with
  % the exact screen resolution
  
  %% Spatial definitions (screen-size defined)
  scaleFactor = 1e-3; % um to mm. Constant factor
  % halfSizeX,Y: half physical size of the SLM's active area. Taken with 
  % the datasheet parameters
  halfSizeX = monitorSize(1)*pixSize*scaleFactor/2;
  halfSizeY = monitorSize(2)*pixSize*scaleFactor/2;
  % x,y: vectors of SLM's physical size  
  x = linspace(-halfSizeX,halfSizeX,monitorSize(1));                                     
  y = linspace(-halfSizeY,halfSizeY,monitorSize(2)); 
  X = Xcoord2; Y = Ycoord2; % Screen coordinates
  
end

%% Aspect ratio application for the scaling
  % Regarding the drawing of the masks on the SLM screens:
  % Original X (output of f_MakeScreenCoords): the circular mask is drawn 
  % as an ellipse due to the screen transformation (elliptical scaling of 
  % the space)
  % Xrescaled: the spatial scaling is compensated and the circular mask is
  % drawn normally on the whole screen
%   a = 0; % Dummy variable: flag used later
  if circularMask == 1 && plotMask == 2 && ...
     (MaxMask == 1 || coordType == 2) && (maskSel ~= 5 && maskSel ~= 6)
     Xrescaled = AspectRatio*X; % Used for the mask generation: X scaling
     X = Xrescaled; % Circular truncation in full screen
%      a = 1; % Dummy variable: flag used later
  end % Otherwise X=X and one has the elliptical truncation in full screen

%% Polar coordinates for the PC
Xpc = X; Ypc = Y; % Needed for the PC coordinates later on
[phiPC,rPC] = cart2pol(Xpc,Ypc); % Without shifts and no scaling: mask is 
                                 % always circular and centered
                                 
%% Half support for the shift scalling and the mask truncation radius
HalfSupportX = f_MatrixHalfSupport(X);
HalfSupportY = f_MatrixHalfSupport(Y);

%% Mask truncation radius
% This is the radius for the circular (or elliptical) pupil truncation
% The selected minimum takes into account both possible screen 
% configurations: landscaped and portrait
rSize = min(HalfSupportX,HalfSupportY);  
% This is the maximum radius that allows to circumscribe a circle inside 
% a square or inside a rectangle
          
%% Shift of the mask (for the SLM)
switch shiftBool 
 case 0
  shiftX = 0; shiftY = 0; % Shift deactivated   
      
 case 1
  % Test if the center of the mask is inside the truncation
  if max(shiftCart) > 100
    error(['The center of the mask is out of the boundaries of the ' ...
           'image, please select a smaller value for "shiftCart" ']);
  end
  shiftCart = shiftCart/100; % Percentage w.r.t the half size 
  % old shift in cm: shiftCart = spaceSupport*shiftCart/100
  shiftX = shiftCart(2); % Cartesian shift in x
  shiftY = shiftCart(1); % Cartesian shift in y
 
 case 2 % Self-centering algorithm
  % Pending
end

%% Shift scalling so that it is a percentage
% Takes into account the half support of X and Y
shiftX = shiftX*HalfSupportX;
shiftY = shiftY*HalfSupportY;

%% Shift application
% The signs of the shifts account for the cartesian coordinates convention
Xslm = X - shiftX; % Shifted X for the SLM 
Yslm = Y + shiftY; % Shifted Y for the SLM.

%% Polar coordinates for the SLM
% X,Y variables redefined for being used in the EGV and Fork masks

[phiSLM,rSLM] = cart2pol(Xslm,Yslm); % Polar coordinates with an added
                                     % shift. The signs compensate the 
                                     % normal cartesian convention for 
                                     % displacing the phase mask
                             
%% Zernike polynomials square size
ZernikeSize = min(monitorSize); % Minumum size of the screen for the
                                % square-size Zernike polynomials
 
end