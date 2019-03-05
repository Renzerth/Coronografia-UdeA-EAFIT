%% Generate a Zernike phase mask (wavefront)
function [mask,wrapMask,wrapMaskFig] = f_ZernikeMask(X,Y,r,z_coeff,z_a,z_frac, ...
L,gl,phaseValues,mingl,maxgl,levShft,z_pupil,z_disp_wrap,z_plot,normMag, ...
binMask,binv,rSize,monitorSize,scrnIdx,coordType,MaxMask,plotMask)
% Characterizes the aberrations of the system
% Inputs:
%  X,Y: A grid of the spatial vector: 2D Cartesian coordiantes. They
%  may already have the shiftCart
%  r: polar coordinate
%  z_coeff: 
%     selected aberrations: any combination on [1,20]; zernike
%     coefficients. ANSI standard. Zernike with desired weights.
%     Put the aberration j-numbers: j = [0 1 2 3 4 ... 20]
%     A) Random surface: -1; Zernike with random weights on [1,20]
%     B) Column vector with normalized coefficients on [0,1]
%     C) Row vector with the j #s that one wants to plot (only integers)
%
%     Examples:
%     A) z_coeff = -1 (just that)
%     B) z_coeff = [0.1 -0.13 0.15 0 -1 0.78]' (transposed)
%     C) z_coeff = [2 3 1 5 7] (normal)
%  z_a: arbitrary constant; the bigger, the more intense the mask;ref: a=20
%  z_frac: to adjust the wrapped phase; ref: 0.125
%  L: laser wavelength [um]
%  gl: gray levels of Zernike. Normally 256
%  phaseValues: discretized phi vector on [-pi,pi].
%  mingl,maxgl: minimum/maximum gray level depth. Ref: 0,255
%  levShft: corresponds to the brightness or constant shift of the gl's
%  z_pupil: defines pupil relative size (w.r.t. sSize), like a percentage
%  z_disp_wrap: unwrapped (0) or wrapped mask (1)
%  z_plot: plot (1); no plot (0)
%  normMag: normalize magnitude. yes(1); no(0)
%  binMask: binarizes the mask w.r.t the max and min of the phase (boolean)
%  binv: binary inversion of the mask: yes(1); no(0). Only applies when 
%        binMask=1. It is usefull to be applied for odd p's on LG beams
%  rSize: radius for the circular (or elliptical) pupil truncation
%  monitorSize: size of the selected screen for coordType = 2 or of the 
%  grid (sSize) for coordType = 1  
%  scrnIdx: screen number selector. In [1,N] with N the # of screen
%  coordType: type of calculation of the spatial coordinates. def: 2 
%    -1: size defined by the user, space support defined by the SLM to use
%    -2: size defined by the resolution of the selected screen    
%  MaxMask: defines if the mask should be maximized when coordType = 1
%           -0: custom-size mask that depends on the variable sSize   
%           -1: maximizes the mask for coordType = 1
%           -2: maximized mask but keeping its rectangular fashion
%  plotMask:  no (0); on the screen (1); on the SLM (2); on the screen, but
%             a surface (3)
%
% Outputs:
%  Zernike phase mask. Complex structure that has not been truncated and is
%  wrapped on [-pi,pi]. mask = exp(i*UnwrappedMask).
%  wrapMask: truncation and angle operations on mask.
%  wrapMaskFig: figure handler if needed outside the function
%
% Missing:
%  it still doesn't modulate on 2pi completely

%% Parameters (input z_coeff vector verification)
o1 = size(z_coeff); % Dimensions of z_coeff
o2 = abs(z_coeff) > ones(size(z_coeff)); % Ask if elemts are bigger than one
o2 = sum(o2); % Sum all elements; zero if they are all smaller than one
if z_coeff(1) == -1 && isscalar(z_coeff) == 1 % A); when z_coeff = -1
 z_vec = 2*pi*z_frac*z_a/L*randn(1,15)'; % Transposed random vector of integers
elseif o1(1) == 1 &&  sum(mod(z_coeff,1))== 0 % [#s] and integers % C)
 z_vec = zeros(1,15)'; % Transposed vector
 z_coeff =  z_coeff + 1; % Correction in order to get correctly the
                         % polynomials (j-th)
 z_vec(z_coeff) = 2*pi*z_frac*z_a/L; % Zernike weight vect: characterizes
                                 % the aberrations to be plotted
elseif o1(1) == 1 &&  o2 == 0 % Column vector and #s less than 1 [#s]' % B)
 z_vec = zeros(1,15)'; % Transposed vector
 for i = 1:size(z_coeff,2)
      z_vec(i) = 2*pi*z_frac*z_a*(z_coeff(i))/L;
 end
else
  warning('Not valid input');
  z_vec = zeros(1,2)'; % Initialization so that there are no errors
end

%% Zernike Polynomials
unwrapmask = f_ZernikeBuilder(X,Y,z_vec,z_pupil,z_plot); 
% (vector, pupil size, Matrix size (zernike phase size), graph:1 or not:0)
mask = exp(1i*unwrapmask); % Wrapped mask

%% Plot the mask
if z_disp_wrap == 1 % Wrapped phase
    abs_ang = 2; % Phase plot
    tit ='Wrapped phase mask';
    % mask = mask; % No change
    str = ''; % Empty, it only works for abs_ang = 0
else % disp_wrap = 0
    abs_ang = 0; % Custom input in order to not wrap the phase
    tit = 'Unwrapped phase mask'; % The amplitude title is replaced  
    mask = unwrapmask; % Unwrapped mask
    str = 'Unwrapped phase value';
end
MaskPupil = 0; % Always, since Zernike has its own pupil
[wrapMask,wrapMaskFig] = f_ProjectMask(r,mask,gl,phaseValues,mingl,...
maxgl,levShft,normMag,binMask,binv,MaskPupil,rSize,monitorSize,scrnIdx,tit, ...
str,coordType,abs_ang,MaxMask,plotMask);

%% Phase wrapping (optional test for curiosity) 
% The wrapped phase should do complete cycles of 2pi
% min(b,[],'omitnan') % Negative
% min(cos(b(:)),[],'omitnan') % positive since the phase was wrapped 
                              % for intensity   
                              
end                              