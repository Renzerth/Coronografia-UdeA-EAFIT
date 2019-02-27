function [x,y,Xslm,Yslm,rSLM,phiSLM,Xpc,Ypc,rPC,phiPC,ZernikeSize, ...
monitorSize] = f_DefineSpace(sSupport,sSize,shiftCart,shiftBool,pixSize,...
scrnIdx,circularMask,coordType,MaxMask,plotMask,maskSel)
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
%  circularMask: Only works when coordType = 2 and when maskSel ~= (5,6)
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
%
% Outputs:
% x,y: cartesian coordinates
% Xslm,Yslm,rSLM,phiSLM,Xpc,Ypc,rPC,phiPC: meshgrids of 2D cartesian and 
%                                          polar coordinates of both the 
%                                          SLM and the PC
% ZernikeSize: screen size for the Zernike polynomials generation
% monitorSize: resolution of the monitor selected by the scrnIdx

%% Screen coordinates (for coordType=2) and Aspect Ratio for a general use
  enablechange = false; 
  % false: won't change default figure display monitor. Leave this value as
  % zero always as the figure display monitor will be changed later on.
  % This is changed if one wants to display in the SLMs when plotMask = 2
  [Xcoord2,Ycoord2,AspectRatio,monitorSize] = ...
  f_MakeScreenCoords(scrnIdx,enablechange); % Calculates the monitor size

%% Coordinate type selection  
switch coordType
 case 1 % Size defined by the user, space support defined by the SLM to use
  % This applies when one won't project a full screen mask and a desired 
  % image resolution is wanted 
  %% Spatial definitions (size is user defined)
  sSupport = sSupport/2; % Half support of the SLM window in cm
  spaceVector = -sSupport:2*sSupport/(sSize-1):sSupport; 
  % Symmetric space
  [Xcoord1,Ycoord1] = meshgrid(spaceVector); % A symmetric grid: 2D Cartesian 
                                 % coordinates
  x = spaceVector; % Cartesian x-vector
  y = x; % Cartesian y-vector: square grid
  % Xslm = X;
  monitorSize = size(Xcoord1); % Square coordinates: [sSize,sSize]
  
  enablechange = false; % Always false
  [~,~,AspectRatio,~] = f_MakeScreenCoords(scrnIdx,enablechange);
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
  Xrescaled = AspectRatio*X; % Used for the mask generation: scaling of X
  if circularMask == 1 && plotMask == 2 && ...
     (MaxMask == 1 || coordType == 2) && (maskSel ~= 5 && maskSel ~= 6)
     X = Xrescaled; % Circular truncation in full screen
  end % Otherwise X=X and one has the elliptical truncation in full screen

%% Polar coordinates for the PC
Xpc = X; Ypc = Y; % Needed for the PC coordinates later on
[phiPC,rPC] = cart2pol(Xpc,Ypc); % Without shifts and no scaling: mask is 
                                 % always circular and centered
                          
%% Shift of the mask (for the SLM)
switch shiftBool 
 case 0
  shiftX = 0; shiftY = 0; % Shift deactivated   
      
 case 1
  shiftCart = shiftCart/100; % Percentage w.r.t the half size 
  % old shift in cm: shiftCart = spaceSupport*shiftCart/100
  shiftX = shiftCart(2); % Cartesian shift in x
  shiftY = shiftCart(1); % Cartesian shift in y
 
 case 2 % Self-centering algorithm
  % Pending
end

%% Polar coordinates for the SLM
% X,Y variables redefined for being used in the EGV and Fork masks
% The signs of the shifts account for the cartesian coordinates convention
Xslm = X - shiftX; % Shifted X for the SLM
Yslm = Y + shiftY; % Shifted Y for the SLM.
[phiSLM,rSLM] = cart2pol(Xslm,Yslm); % Polar coordinates with an added
                                     % shift. The signs compensate the 
                                     % normal cartesian convention for 
                                     % displacing the phase mask
                             
%% Zernike polynomials square size
ZernikeSize = min(monitorSize); % Minumum size of the screen for the
                                % square-size Zernike polynomials
end