function [rSize,x,y,Xslm,Yslm,rSLM,phiSLM,Xpc,Ypc,rPC,phiPC,monitorSize,...
shiftCart,mainLyotRadius] = f_DefineSpace(vidSelfCent,sSupport,sSize, ...
shiftCart,shiftMask,PP,pixSize,apRad,scrnIdx,circularMask,z_pupil, ...
coordType,MaxMask,SLMcenterWisdom,camera,cameraPlane,exposure, ...
format,fps,measSimulated,maskSel)
% Inputs:
%  vid: video input object
%  sSupport: full side-length of the SLMs (or unitary without an SLM)
%  sSize: bits for grey levels; 2^k is the resolution (size of x and y)
%     Default: 10. Size is calculated as 2^k - 1
%     Only works when coordType = 1
%  shiftCart: [xshift,yshift], works when shiftMask = 1. Percentages of
%             movement of the HALF size of the screen (cartesian coordinates
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
%             shiftMask=0 or with the camera-compensated shift when 
%             shiftMask=2 (units for all cases: percentage in [0,1], 
%             relative)
%
% Outputs:
%  rSize: radius for the circular (or elliptical) pupil truncation
%  x,y: cartesian coordinates
%  Xslm,Yslm,rSLM,phiSLM,Xpc,Ypc,rPC,phiPC: meshgrids of 2D cartesian and 
%                                           polar coordinates of both the 
%                                           SLM and the PC
%  slm: shifts possibly applied; pc: shifts are never applied
%  ZernikeSize: screen size for the Zernike polynomials generation
%  monitorSize: resolution of the monitor selected by the scrnIdx

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
%      Xrescaled = AspectRatio*X; % Used for the mask generation: X scaling
%      X = Xrescaled; % Circular truncation in full screen
     Yrescaled = Y/AspectRatio;
     Y = Yrescaled;
     flagAR = 1; % Flag for the AspectRatio application
  end % Otherwise X=X and one has the elliptical truncation in full screen
                                 
%% Half support for the shift scaling and the mask truncation radius
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
widomexist = (exist(SLMcenterWisdom,'file') == 2);
if isempty(vidSelfCent) && shiftMask == 2 && ~(widomexist)
  error(['Set measSimulated = 0 and meas = 1 in order to use' ...
    ' shiftMask = 2. shiftMask should be set to zero or one.']);
elseif isempty(vidSelfCent) && shiftMask == 2  && measSimulated == 1
  warning('SLMwisdom.mat will be loaded for measSimulated = 1')  
end

switch shiftMask 
 case 0
  shiftX = 0; shiftY = 0; % Shift deactivated   
  shiftCart = [shiftX shiftY];
  
  mainLyotRadius = nan;
  
 case 1 % User-given
  % Test if the center of the mask is inside the truncation
  if max(shiftCart) > 100
    error(['The center of the mask is out of the boundaries of the ' ...
           'image, please select a smaller value for "shiftCart" ']);
  end
  shiftCart = shiftCart/100; % Percentage w.r.t the half size 
  shiftX = shiftCart(1); % Cartesian shift in x
  shiftY = shiftCart(2); % Cartesian shift in y
  shiftCart = [shiftX*AspectRatio, shiftY];
  
  mainLyotRadius = nan;
  
 case 2 % Self-centering algorithm
  shiftCartfine = shiftCart; % User-given for a fine adjustment
  shiftCartfine = shiftCartfine/100; % Percentage w.r.t the half size 
  shiftXfine = shiftCartfine(2)*AspectRatio; % Cartesian shift in x
  shiftYfine = shiftCartfine(1); % Cartesian shift in y
  % The camera here must be Lyot for the 2019's setup
  if ~(exist(SLMcenterWisdom,'file') == 2) % The self centering data
                                           % doesn't exist in Data (folder)
                                           
    %%% Self centering parameters
    PP = PP*1e-6;  % um to m
    lensDiameter = 2*apRad*1e-2; % cm to m
    
    %%% Preview parameters for the self centering
    [N,D] = rat(exposure);
    tit = strcat('Preview of the',{' '},cameraPlane,{' '}, ...
    'camera (',camera,')',{' '},'[Exposure:',{' '},num2str(N),'/', ...
    num2str(D),';',{' '},'format:',{' '},num2str(format),';',{' '}, ...
    'fps:',{' '},num2str(fps),']');
    
   %% Self-centering algorithm
    % these outputs are not used in the meantime (2019): 
    % [shiftY,shiftX,systemPupilPixelSize,mainDataCenter,mainDataRadius] 
    [shiftX, shiftY,~,~,mainLyotRadius] = ...
    f_selfCenterSLMmask(PP,lensDiameter,scrnIdx,vidSelfCent,tit);
    if coordType == 2
      [shiftX,shiftY] = f_calcHScoorToSgnCoor(shiftX/monitorSize(1), ...
        shiftY/monitorSize(2));
    end
    shiftCart = [shiftX, shiftY];
    
    %% Save SLMwisdom.mat
    % Explanation: save(directory+filename,variables)
    save(SLMcenterWisdom,'shiftCart','mainLyotRadius')
    
  else % The self centering data is available in Data (folder)
    %% Load SLMwisdom.mat
    load(SLMcenterWisdom,'shiftCart','mainLyotRadius');
    shiftX = shiftCart(1); shiftY = shiftCart(2);
    
    %% Fix data format incompatibility
    if coordType == 2 && max(abs(shiftCart)) >= 1
      % This converts the shifts from coordType 1 to 2
      [shiftX,shiftY] = f_calcHScoorToSgnCoor(shiftX/monitorSize(1), ...
        shiftY/monitorSize(2));
    elseif coordType == 1 && max(abs(shiftCart)) <= 1
      % This converts the shifts from coordType 2 to 1
      [shiftX,shiftY] = f_calcSgnCoorToHScoor(shiftX, shiftY);
      % OLD:
      %       shiftY = round(shiftY*monitorSize(2));
      %       shiftX = round(shiftX*monitorSize(1));
    end
    shiftCart = [shiftX, shiftY];
  end
 end % switch shiftMask
 % After this switch, shiftCart is taken as the output so that all the
 % masks 
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
if shiftMask ~= 2 
  shiftX = shiftX*HalfSupportX; % Shift converted to the HalfSupport's units
  shiftY = shiftY*HalfSupportY; % Shift converted to the HalfSupport's units
end

%% Shift application

if flagAR == 1 % Aspect ratio already applied on Y
  shiftY = shiftY/AspectRatio;
end

if shiftMask == 2 
    % The signs of the shiftX,Y account for "perhapsAworkingDEMO.m"
    % convention and the signs of the shiftX,Yfine account for the 
    % cartesian coordinate convention (as for shiftMask = [0,1])
    Xslm = X - shiftX + shiftXfine*shiftX; % shiftXfine*shiftX: a percent.
                                           % of the current shift
    Yslm = Y - shiftY - shiftYfine*shiftY;
else % shiftMask = [0,1]
    % The signs of the shifts account for the cartesian coordinate
    % convention
    Xslm = X + shiftX;
    Yslm = Y - shiftY;
end

%% Polar coordinates for the SLM
% X,Y variables redefined for being used in the EGV and Fork masks
[phiSLM,rSLM] = cart2pol(Xslm,Yslm); % Polar coordinates with an added
                                     % shift. The signs compensate the 
                                     % normal cartesian convention for 
                                     % displacing the phase mask
% [Xr,Yr,~,~] = f_MakeScreenCoords(3,false);
% Xslm = Xr - shiftX;
% Yslm = (Yr - shiftY)/AspectRatio;
% [phiSLM,rSLM] = cart2pol(Xslm, Yslm);                                     
%% Termination of the hardware 
f_releaseCamera(vidSelfCent); % Clear vid for the self-centering stage

end