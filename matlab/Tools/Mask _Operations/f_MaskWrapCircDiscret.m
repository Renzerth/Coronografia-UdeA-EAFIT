function [wrapMask,customMap] = f_MaskWrapCircDiscret(r,mask, ...
phaseValues,binMask,binv,MaskPupil,rSize,plotMask)
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
%  rSize: radius for the circular (or elliptical) pupil truncation
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

[wrapMask,customMap] = f_discretizeMask(phaseValues,wrapMask); % Mask 
                                                           % discretization

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
  
  %% Mask padarray with zeros if needed
  % Only used for Zernike masks (the only one assumed to be generated with 
  % a squared-size):
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

