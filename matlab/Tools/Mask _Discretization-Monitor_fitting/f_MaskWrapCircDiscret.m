function wrappedMask = f_MaskWrapCircDiscret(r,mask,binMask,binv,glphi,...
                       mingl,maxgl,levShft,coordType,plotMask)
% Multiplies the phase mask by the maximum circle size with its outer
% borders containing the minimum value of the phase (normally -pi)
% Wraps the phase with the function "angle"
% Discretizes the mask with specific gray levels and dynamic range
%
% Inputs:
%  r: polar coordinate (in cm)
%  mask: complex structure that has not been truncated and is wrapped on
%        [-pi,pi]. mask = exp(i*UnwrappedMask) 
%  binMask: binarizes the mask w.r.t the max and min of the phase (boolean)
%  binv: binary inversion of the mask: yes(1); no(0). Only applies when 
%        binMask=1. It is usefull to be applied for odd p's on LG beams
%  glphi: discretized phi vector on [-pi,pi].
%  mingl,maxgl: minimum/maximum gray level depth. Ref: 0,255
%  levShft: corresponds to the brightness or constant shift of the gl's
%  coordType: type of calculation of the spatial coordinates. def: 2 
%    -1: size defined by the user, space support defined by the SLM to use
%    -2: size defined by the resolution of the selected screen    
%  plotMask:  no (0); on the screen (1); on the SLM (2); on the screen, but
%             a surface (3)
%
% Output:
%  wrappedMask: real-valued mask (angle), wrapped on [-pi,pi] and truncated

%% Mask wrapping
wrappedMask = angle(mask); % Phase of the mask on [-pi, pi]. Real-valued

%% Discretized phase mask
% Important: this is applied after the exp(i*mask) was created and then it
% was wrapped with angle so that it is discretized on [-pi,pi]

wrappedMask = f_discretizeMask(wrappedMask,glphi); % Mask discretization
wrappedMask = f_ScaleMatrixData(wrappedMask,mingl,maxgl) + levShft; 
% Scaling to uint8 values

%% Mask Binarization and binary inversion
% binarizes the mask w.r.t the max,mid and min of the phase (boolean)
minMask = min(wrappedMask(:)); % Minimum value of the mask
maxMask = max(wrappedMask(:)); % Maximum value of the mask
midMask = (minMask+maxMask)/2; % Medium value of the mask

if binMask == 1 
    wrappedMask = double(wrappedMask < midMask); 
    % Binarizes on [0,1] by the phase above and bellow the mid value
    % wrappedMask = wrappedMask*minMask; % Establishes the lower binary
    % value as the lowest value of the phase
    if binv == 1 % Binary inversion
        wrappedMask = imcomplement(wrappedMask); % Applied to the angle
    end
end

% In both cases of the next if-else, rMax is found as the maximum radius
% that allows to circumscribe a circle inside an square for coordType == 1
% or inside a rectangle for coordType == 2
if coordType == 1 % User-defined
    rMax = max(r(:));  % the maximum value of r (diagonal of the square)
    rSize = rMax/sqrt(2); % Equals this since twice rSize^2 equals
                          % rmax^2 (Pythagorean theorem)
else % coordType == 2 % Screen-resolution defined
    rSizeVect = size(r); % 2D vector with the size of r
    rMidVect = floor((rSizeVect+1)/2); % mid points of the size of r
    rMove = max(rMidVect); % Maximum midpoint in order to move here
    idx = find(rMidVect == rMove); % finds the index or rMove to determine
                                   % if one should move in the x or the y 
                                   % direction
    if idx == 1 % (for landscape monitors)
        rMax = r(rMove,1); % The rMax is in the y direction 
    else % idx == 2 % (for portrait monitors)
        rMax = r(1,rMove); % The rMax is in the x direction
    end
    % if one has a unitary space, rMax = 1 always
    rSize = rMax; % Both rmax and rsize are equal
end

%% Mask padarray with zeros if needed
% Only used for Zernike masks (the only one assumed to be generated with a
% square-size):
A = size(r) - size(wrappedMask); % Size comparison
if any(A) % True when tests whether any of the elements along various
          % dimensions of an array are nonzero
    if  plotMask == 2 % SLM
    
        idx = find(A ~= 0);
        pad = A(idx);
        if isscalar(pad)
            %% Orientation of the padding
            % Horizontal Padding:
            a1 = [0 pad/2];
            a2 = [0 (pad+1)/2];
            a3 = [0 (pad+1)/2 - 1];
            if idx == 1 % Vertical Padding: same as horizontal but shifted
               a1 = circshift(a1,1); % Circular shift
               a2 = circshift(a2,1); % Circular shift
               a3 = circshift(a3,1); % Circular shift
            end 
            %% Padding
            % Replicated since that the outer value of the Zernike are
            % constant
            if mod(pad,2) == 0 % Even pad (symmetric)
                wrappedMask = padarray(wrappedMask,a1,'replicate','both');
                % A = [zeros(1,pad/2) A zeros(1,pad/2)];
            else % Odd pad (asymmetric)
                wrappedMask = padarray(wrappedMask,a2,'replicate','pre');
                wrappedMask = padarray(wrappedMask,a3,'replicate','pos');
                % A = [zeros(1,(pad+1)/2) A zeros(1,(pad+1)/2-1)]; 
            end  
        else % pad is a vector
            % MISSING
        end
    else % plotMask == 0 or 1 or 3 % PC
        % "Anti" pad array: r takes the mask's size
        mYmid = floor((size(wrappedMask,1)+1)/2);
        mXmid = floor((size(wrappedMask,2)+1)/2);
        rXmid = floor((size(r,2)+1)/2);
        rYmid = floor((size(r,1)+1)/2);
        r = r(rYmid-mYmid+1:rYmid+mYmid,rXmid-mXmid+1:rXmid+mXmid); 
    end
end

%% Phase mask times a circular (or elliptical) aperture
% Depends on the variable circularMask inside f_DefineSpace.m
binCirc = double(r <= rSize); % Binary mask.
wrappedMask = wrappedMask.*binCirc; % Binary mask. Range: [minMask,maxMask]
wrappedMask(r > rSize) = nan; % Outside the circular (elliptical) mask
% OLD: wrappedMask(r > rSize) = minMask; % Outside the circular pupil puts
% the smallest value of the phase mask

end

