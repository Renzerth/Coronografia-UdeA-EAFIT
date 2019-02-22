switch coordType
 case 1 % Size defined by the user, space support defined by the SLM to use
  %% Spatial definitions
  sSize = 2^k-1; % Number of samples; odd number so that vortex gets
                 % centered (spatial size); Spatial size. ref: 2^k-1
  spaceSupport = spaceSupport/2; % Half support of the SLM window in cm
  spaceVector = -spaceSupport:2*spaceSupport/(sSize-1):spaceSupport; 
  % Symmetric space
  [X,Y] = meshgrid(spaceVector); % A symmetric grid: 2D Cartesian 
                                 % coordinates
  x = spaceVector; % Cartesian x-vector
  y = x; % Cartesian y-vector: square grid

 case 2 % Size defined by the resolution of the selected screen                 
  %% Screen coordinates
  enablechange = true; 
  % false: won't change default figure display monitor. Leave this value as
  % zero always as the figure display monitor will be changed later on.
  % This is changed if one wants to display in the SLMs when plotMask = 2
  [X,Y,aspectRatio,monitorSize] = ...
                             f_makeScreenCoordinates(scrnIdx,enablechange);
  % Calculates the monitor size
  scaleFactor = 1e-3; % um to mm. Constant factor
  % halfSizeX,Y: half physical size of the SLM's active area. Taken with 
  % the datasheet parameters
  halfSizeX = monitorSize(1)*pixSize*scaleFactor/2;
  halfSizeY = monitorSize(2)*pixSize*scaleFactor/2;
  % x,y: vectors of SLM's physical size  
  x = linspace(-halfSizeX,halfSizeX,monitorSize(1));                                     
  y = linspace(-halfSizeY,halfSizeY,monitorSize(2)); 
end

%% Polar coordinates with a shift of the mask
if shiftBool == 1
 shiftCart = shiftCart/100; % Percentage w.r.t the half size 
 % old shift in cm shiftCart = spaceSupport*shiftCart/100
 shiftX = shiftCart(2); % Cartesian shift in x
 shiftY = shiftCart(1); % Cartesian shift in y
else
 shiftX = 0; shiftY = 0; % Shift deactivated   
end
[phi,r] = cart2pol(X-shiftX,(Y+shiftY)); % Polar coordinates with an added
                                         % shift. The signs compensate the 
                                         % normal cartesian convention for 
                                         % displacing the phase mask