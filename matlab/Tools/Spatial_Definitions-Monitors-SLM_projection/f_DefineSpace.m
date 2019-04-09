function [rSize,x,y,Xslm,Yslm,rSLM,phiSLM,Xpc,Ypc,rPC,phiPC,monitorSize,...
shiftCart] = f_DefineSpace(vid,sSupport,sSize,shiftCart,shiftMask,PP, ...
pixSize,apRad,scrnIdx,circularMask,z_pupil,coordType,MaxMask, ...
SLMcenterWisdom,maskSel)
% Inputs:
%  vid: video input object
%  sSupport: full side-length of the SLMs (or unitary without an SLM)
%  sSize: bits for grey levels; 2^k is the resolution (size of x and y)
%     Default: 10. Size is calculated as 2^k - 1
%     Only works when coordType = 1
%  shiftCart: [yshift,xshift], works when shiftMask = 1. Percentages of
%             movement of the total size of the mask (cartesian coordinates
%             convention). Ranges per shift: [0,100] (percentage)  
%  shiftMask: only shifts when plotMask = 2
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
%  z_pupil: defines pupil relative size (w.r.t. sSize), like a percentage
%  coordType: type of calculation of the spatial coordinates. def: 2 
%    -1: size defined by the user, space support defined by the SLM to use
%    -2: size defined by the resolution of the selected screen
%  MaxMask: defines if the mask should be maximized when coordType = 1
%           -0: custom-size mask that depends on the variable sSize   
%           -1: maximizes the mask for coordType = 1
%           -2: maximized mask but keeping its rectangular fashion
%  SLMcenterWisdom: directory to store shiftCart
%  maskSel: selects a specific mask
%  shiftCart: same as input but over 100 and filled with zeros when
%                 shiftMask=0 or with the camera-compensated shift when 
%                 shiftMask=2 (units for all cases: percentage in [0,1])
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

%% Coordinate type selection  
switch coordType
 case 1 % Size defined by the user, space support defined by the SLM to use
  % This applies when one won't project a full screen mask and a desired 
  % image resolution is wanted 
  
  %% Spatial definitions (size is user defined)
  % Physical size except for VPL and Zernike, as they have to be in a 
  % unitary space
  if ismember(maskSel,[2,5,6])
    sSupport = 2; % Unitary space since one makes sSupport/2 = 1 here
  end
  sSupport = sSupport/2; % Half support of the SLM window in cm
  spaceVector = -sSupport:2*sSupport/(sSize-1):sSupport; 
  % Symmetric space
  [Xcoord1,Ycoord1] = meshgrid(spaceVector); % A symmetric grid: 2D 
                                             % Cartesian coordinates
  x = spaceVector; % Cartesian x-vector
  y = x; % Cartesian y-vector: square grid
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
  % X,Y are not created with x,y. X,Y are on [-1,1]. x,y: physical size
  X = Xcoord2; Y = Ycoord2; % Screen coordinates
  
end

%% Polar coordinates for the PC (plotMask=1,3)
Xpc = X; Ypc = Y; % Needed for the PC coordinates later on
[phiPC,rPC] = cart2pol(Xpc,Ypc); % Without shifts and no scaling: mask is 
                                 % always circular and centered

%% Aspect ratio application for the scaling for the SLM (plotMask=2)
  % Regarding the drawing of the masks on the SLM screens:
  % Original X (output of f_MakeScreenCoords): the circular mask is drawn 
  % as an ellipse due to the screen transformation (elliptical scaling of 
  % the space)
  % Xrescaled: the spatial scaling is compensated and the circular mask is
  % drawn normally on the whole screen
  flagAR = 0; % Flag for the AspectRatio application (not applied)
  if circularMask == 1 && (MaxMask == 1 || coordType == 2) 
     Xrescaled = AspectRatio*X; % Used for the mask generation: X scaling
     X = Xrescaled; % Circular truncation in full screen
     flagAR = 1; % Flag for the AspectRatio application
  end % Otherwise X=X and one has the elliptical truncation in full screen
                                 
%% Half support for the shift scalling and the mask truncation radius
HalfSupportX = f_MatrixHalfSupport(X);
HalfSupportY = f_MatrixHalfSupport(Y);
% Both are 1 for coordType = 2

%% Mask truncation radius
% This is the radius for the circular (or elliptical) pupil truncation
% The selected minimum takes into account both possible screen 
% configurations: landscape and portrait
rSize = min(HalfSupportX,HalfSupportY);  
% This is the maximum radius that allows to circumscribe a circle inside 
% a square or inside a rectangle
          
%% Shift of the mask (for the SLM: plotMask=2)
% Check first if vid exists, otherwise set shiftMask to zero
if isempty(vid) && shiftMask == 2
  error(['Set measSimulated = 0 and meas = 1 in order to use' ...
    ' shiftMask = 2. shiftMask should be set to zero.']);
end

switch shiftMask 
 case 0
  shiftX = 0; shiftY = 0; % Shift deactivated   
  shiftCart = [shiftX shiftY];    
  
 case 1 % User-give
  % Test if the center of the mask is inside the truncation
  if max(shiftCart) > 100
    error(['The center of the mask is out of the boundaries of the ' ...
           'image, please select a smaller value for "shiftCart" ']);
  end
  shiftCart = shiftCart/100; % Percentage w.r.t the half size 
  shiftX = shiftCart(2); % Cartesian shift in x
  shiftY = shiftCart(1); % Cartesian shift in y
  shiftCart = [shiftY shiftX];
  
 case 2 % Self-centering algorithm

  %% Self Centering algorithm
  % The camera here must be Lyot for the 2019's setup
  if ~(exist(SLMcenterWisdom,'file') == 2) % The self centering data
                                           % doesn't exist in Data (folder)
    PP = PP*1e-6;  % um to m
    lensDiameter = 2*apRad*1e-2; % cm to m
    [shiftY, shiftX,~,~,~] = ...
    f_selfCenterSLMmask(PP, lensDiameter, scrnIdx, vid);
    [shiftX,shiftY] = f_calcHScoorToSgnCoor(shiftX/monitorSize(1), ...
                                            shiftY/monitorSize(2));
    shiftCart = [shiftY, shiftX];
    % Explanation: save(directory+filename,variables)
    save(SLMcenterWisdom,'shiftCart')
  else % The self centering data is available in Data (folder)
    load(SLMcenterWisdom,'shiftCart');
  end
end
% After this switch, shiftCart is taken as the output so that all the masks
% shown during the measurement have this shift

%% Check the shift for the Zernike polynomials
% Right now, shiftX and shiftY are percentages but z_pupil scaled by the
% half supports if it is applied above (so that flagAR = 1)
if flagAR == 1
  z_pupilX = z_pupil/HalfSupportX; % Half support scaling in X
  z_pupilY = z_pupil/HalfSupportY; % Half support scaling in Y
else
   z_pupilX = z_pupil; % z_pupil is not scaled in X
   z_pupilY = z_pupil; % z_pupil is not scaled in Y
end 

cond1 = (abs(shiftX) + z_pupilX) > 1; % X shift percentage constraint
cond2 = (abs(shiftY) + z_pupilY) > 1; % Y shift percentage constraint
if (cond1 || cond2) && (maskSel == 5 || maskSel == 6) 
  % Error for the SLM coordinates and Zernike masks
  % Zernike masks must always be circular and never truncated
  disp(['For the selected z_pupil, shiftX should be at maximum ' ...
      num2str((1-z_pupilX)*100) ' and shiftY ' num2str((1-z_pupilY)*100)]);
  error(['Shifts out of bounds for the Zernike polynomials: the whole ' ...
        'pupil must lie inside the screen in order for them to be ' ...
        'almost orthogonal. Please set "shiftCart" or "z_pupil" smaller']); 
end
  
%% Shift scaling so that it is a percentage
% Takes into account the half support of X and Y and a percentage of them
% is taken
shiftX = shiftX*HalfSupportX; % Shift converted to the HalfSupport's units
shiftY = shiftY*HalfSupportY; % Shift converted to the HalfSupport's units

%% Shift application
% The signs of the shifts account for the cartesian coordinates convention
Xslm = X - shiftX; % Shifted X for the SLM 
Yslm = Y + shiftY; % Shifted Y for the SLM

%% Polar coordinates for the SLM
% X,Y variables redefined for being used in the EGV and Fork masks
[phiSLM,rSLM] = cart2pol(Xslm,Yslm); % Polar coordinates with an added
                                     % shift. The signs compensate the 
                                     % normal cartesian convention for 
                                     % displacing the phase mask
      
%% Termination of the hardware clear vid for any stage of this algorithm
f_releaseCamera(vid);              

end