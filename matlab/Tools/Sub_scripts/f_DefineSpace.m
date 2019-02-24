function [x,y,Xslm,Yslm,rSLM,phiSLM,Xpc,Ypc,rPC,phiPC,sSize,monitorSize] = ...
f_DefineSpace(spaceSupport,k,shiftCart,pixSize,scrnIdx,circularMask, ...
shiftBool,coordType)
switch coordType
 case 1 % Size defined by the user, space support defined by the SLM to use
  % This applies when one won't project a full screen mask and a desired 
  % image resolution is wanted 
  %% Spatial definitions
  sSize = 2^k - 1; % Number of samples; odd number so that vortex gets
                   % centered (spatial size); Spatial size. ref: 2^k-1
  spaceSupport = spaceSupport/2; % Half support of the SLM window in cm
  spaceVector = -spaceSupport:2*spaceSupport/(sSize-1):spaceSupport; 
  % Symmetric space
  [X,Y] = meshgrid(spaceVector); % A symmetric grid: 2D Cartesian 
                                 % coordinates
  x = spaceVector; % Cartesian x-vector
  y = x; % Cartesian y-vector: square grid
  % Xslm = X;
  monitorSize = size(X);% Square coordinates

 case 2 % Size defined by the resolution of the selected screen     
  % This applies when a full screen mask will be displayed for the SLM with
  % the exact screen resolution
  %% Screen coordinates
  enablechange = false; 
  % false: won't change default figure display monitor. Leave this value as
  % zero always as the figure display monitor will be changed later on.
  % This is changed if one wants to display in the SLMs when plotMask = 2
  [X,Y,AspectRatio,monitorSize] = f_MakeScreenCoords(scrnIdx,enablechange);
  % Calculates the monitor size
  scaleFactor = 1e-3; % um to mm. Constant factor
  % halfSizeX,Y: half physical size of the SLM's active area. Taken with 
  % the datasheet parameters
  halfSizeX = monitorSize(1)*pixSize*scaleFactor/2;
  halfSizeY = monitorSize(2)*pixSize*scaleFactor/2;
  % x,y: vectors of SLM's physical size  
  x = linspace(-halfSizeX,halfSizeX,monitorSize(1));                                     
  y = linspace(-halfSizeY,halfSizeY,monitorSize(2)); 

  %% Aspect ratio application
  % Regarding the drawing of the masks on the SLM screens:
  % Original X (output of f_MakeScreenCoords): the circular mask is drawn 
  % as an ellipse due to the screen transformation (elliptical scaling of 
  % the space)
  % Xrescaled: the spatial scaling is compensated and the circular mask is
  % drawn normally on the whole screen
  Xrescaled = AspectRatio*X; % Used for the mask generation on the pc. It 
                             % is never shifted
  
end

%% Polar coordinates for the PC
Xpc = X; Ypc = Y; % Needed for the PC coordinates later on
[phiPC,rPC] = cart2pol(Xpc,Ypc); % Without shifts and no scaling: mask is 
                                 % always circular and centered
                          
%% Polar coordinates with a shift of the mask (for the SLM)
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

if circularMask == 1 && coordType == 2 % X is used as Xrescaled 
   X = Xrescaled; % Circular truncation in full screen
end % Otherwise X=X and one has the elliptical truncation in full screen

% X,Y variables redefined for being used in the EGV and Fork masks
% The signs of the shifts account for the cartesian coordinates convention
Xslm = X - shiftX; % Shifted X for the SLM
Yslm = Y + shiftY; % Shifted Y for the SLM.
[phiSLM,rSLM] = cart2pol(Xslm,Yslm); % Polar coordinates with an added
                                     % shift. The signs compensate the 
                                     % normal cartesian convention for 
                                     % displacing the phase mask
                             
%% Zernike
sSize = min(min(size(X)),min(size(Y)));
end