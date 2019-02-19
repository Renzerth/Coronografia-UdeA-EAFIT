function wrappedMask = f_mask_circ_angle_gl(r,mask,binMask,glphi,mingl,maxgl,levShft)
% Multiplies the phase mask by the maximum circle size with its outer
% borders containing the minimum value of the phase (normally -pi)
% Wraps the phase with the function "angle"
% Discretizes the mask with specific gray levels and dynamic range
%
% Inputs:
%  r: polar coordinate (in cm)
%  mask: complex phase mask, i.e., mask = exp(1i*unwrapped mask)
%  binMask: binarizes the mask w.r.t the max and min of the phase (boolean)
%  glphi: discretized phi vector on [-pi,pi].
%  mingl,maxgl: minimum/maximum gray level depth. Ref: 0,255
%  levShft: corresponds to the brightness or constant shift of the gl's
%
% Output:
%  wrappedMask: real-valued mask (angle), wrapped on [-pi,pi]

%% Mask wrapping
wrappedMask = angle(mask); % Phase of the mask on [-pi, pi]
minMask = min(wrappedMask(:)); % Just a definition that isn't used

%% Mask Binarization
% binarizes the mask w.r.t the max and min of the phase (boolean)
if binMask == 1 
    wrappedMask = double(wrappedMask < 0); % Binarizes on [0,1] by dividing
                                           % the phase above and bellow 0
    wrappedMask = wrappedMask*minMask; % Establishes the lower binary value 
                                       % as the lowest value of the phase
end

%% Phase mask times a pupil aperture
rmax = max(r(:));  % the maximum value of r (diagonal of the square)
%rSize = rmax/sqrt(2); % Equals this since twice rSize^2 equals
                                 % rmax^2 (Pythagorean theorem)

rSize = rmax;                             
                                 binCirc = double(r <= rSize); % Binary mask.
wrappedMask = wrappedMask.*binCirc; % Binary mask. Range: [minMask,maxMask]
wrappedMask(r > rSize) = nan; % Outside the circular pupil = minMask

end

