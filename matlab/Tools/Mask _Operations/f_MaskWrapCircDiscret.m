function [wrapMask,customMap] = f_MaskWrapCircDiscret(r,mask, ...
phaseValues,binMask,binv,MaskPupil,mingl,maxgl,levShft,coordType,plotMask)
% Multiplies the phase mask by the maximum circle size with its outer
% borders containing the minimum value of the phase (normally -pi)
% Wraps the phase with the function "angle"
% Discretizes the mask with specific gray levels and dynamic range
%
% Inputs:
%  r: polar coordinate (in cm)
%  mask: complex structure that has not been truncated and is wrapped on
%        [-pi,pi]. mask = exp(i*UnwrappedMask) 
%  phaseValues: discretized phi vector on [-pi,pi].
%  binMask: binarizes the mask w.r.t the max and min of the phase (boolean)
%  binv: binary inversion of the mask: yes(1); no(0). Only applies when 
%        binMask=1. It is usefull to be applied for odd p's on LG beams
%  MaskPupil: applies a pupil truncation to the mask: (0): no; (1): yes
%  mingl,maxgl: minimum/maximum gray level depth. Ref: 0,255
%  levShft: corresponds to the brightness or constant shift of the gl's
%  coordType: type of calculation of the spatial coordinates. def: 2 
%    -1: size defined by the user, space support defined by the SLM to use
%    -2: size defined by the resolution of the selected screen    
%  plotMask:  no (0); on the screen (1); on the SLM (2); on the screen, but
%             a surface (3)
%
% Output:
%  wrapMask: real-valued mask (angle), wrapped on [-pi,pi] and truncated

%% Mask wrapping
wrapMask = angle(mask); % Phase of the mask on [-pi, pi]. Real-valued
% wrapMask = mod(mask,2*pi) - pi; % Equivalent operation as angle
% mod(x,a) = x - a*floor(x/a) % Called modfloat
 
%% Discretized phase mask
% Important: this is applied after the exp(i*mask) was created and then it
% was wrapped with angle so that it is discretized on [-pi,pi]

[wrapMask,customMap] = f_discretizeMask(phaseValues,wrapMask);

% wrappedMask = f_discretizeMask(wrappedMask,phaseValues); % Mask 
                                                           % discretization
% wrappedMask = f_ScaleMatrixData(wrappedMask,mingl,maxgl) + levShft; 
% Scaling to uint8 values

%% Mask Binarization and binary inversion
% binarizes the mask w.r.t the max,mid and min of the phase (boolean)
minMask = min(wrapMask(:)); % Minimum value of the mask
maxMask = max(wrapMask(:)); % Maximum value of the mask
midMask = (minMask+maxMask)/2; % Medium value of the mask

if binMask == 1 
    wrapMask = double(wrapMask < midMask); 
    % Binarizes on [0,1] by the phase above and bellow the mid value
    % wrappedMask = wrappedMask*minMask; % Establishes the lower binary
    % value as the lowest value of the phase
    if binv == 1 % Binary inversion
        wrapMask = imcomplement(wrapMask); % Applied to the angle
    end
end

%% Mask truncation
if MaskPupil == 1
  %% Radius for the circular (or elliptical) pupil truncation
  % In both cases of the next if-else, rMax is found as the maximum radius
  % that allows to circumscribe a circle inside an square or inside a 
  % rectangle
  
  % if coordType == 1 % User-defined
  %     rMax = max(r(:));  % the maximum value of r (diagonal of the square)
  %     rSize = rMax/sqrt(2); % Equals this since twice rSize^2 equals
  %                           % rmax^2 (Pythagorean theorem)
  % else % coordType == 2 % Screen-resolution defined
  % end

%   rSizeVect = size(r); % 2D vector with the size of r
%   rMidVect = floor((rSizeVect+1)/2); % mid points of the size of r
%   rMove = max(rMidVect); % Maximum midpoint in order to move here
%   idx = find(rMidVect == rMove); % finds the index or rMove to determine
%                                  % if one should move in the x or the y 
%                                  % direction
%   if isscalar(idx)
%     if idx == 1 % (for landscape monitors)
%           rMax = r(rMove,1); % The rMax is in the y direction 
%     else % idx == 2 % (for portrait monitors)
%           rMax = r(1,rMove); % The rMax is in the x direction
%     end
%   else % idx is a two-row vector, meanning that one has square-sized figures
%     rMax = r(1,rMove)/sqrt(2);
%   end
  % if one has a unitary space, rMax = 1 always
  rSizeVect = size(r);
  rmin = min(rSizeVect);
  rSize =  floor((rmin+1)/2);
%   rSize = rMax; % Both rmax and rsize are equal
  
  %% Mask padarray with zeros if needed
  % Only used for Zernike masks (the only one assumed to be generated with a
  % squared-size):
  % Sizes differ when coordType = 2
  A = size(r) - size(wrapMask); % Size comparison
  if any(A) % True when tests whether any of the elements along various
            % dimensions of an array are nonzero
    if  plotMask == 2 % SLM
      %% wrappedMask matrix padding so that its size fits r
      method = 'replicate'; % For padding: 'replicate', 'symmetric',
                            % 'circular' or a scalar
      wrapMask = f_PadMatrix(wrapMask,r,method);
    else % plotMask == 0 or 1 or 3 % PC
      %% r matrix truncation so that its size fits in wrappedMask
      [mXmid,mYmid] = f_ComputeMatrixMidPoints(wrapMask);
      [rXmid,rYmid] = f_ComputeMatrixMidPoints(r); 
      r = f_TruncateMatrix(wrapMask,mXmid,mYmid,r,rXmid,rYmid);
    end
  end

  %% Phase mask times a circular (or elliptical) aperture
  % Depends on the variable circularMask inside f_DefineSpace.m
  binCirc = double(r <= rSize); % Binary mask.
  wrapMask = wrapMask.*binCirc; % Binary mask. Range: [minMask,maxMask]
  wrapMask(r > rSize) = nan; % Outside the circular (elliptical) mask
  % OLD: wrappedMask(r > rSize) = minMask; % Outside the circular pupil puts
  % the smallest value of the phase mask
end
end

