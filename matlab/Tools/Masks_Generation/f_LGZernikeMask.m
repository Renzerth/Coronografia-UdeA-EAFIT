%% Laguerre Gauss Binary masks + Zernike Phases
function [mask,wrapMask,wrapMaskFig] = f_LGZernikeMask(X,Y,r,phi, ...
phaseValues,tc,s,ph0,p,WsizeRatio,z_coeff,z_a,L,z_pupil, ...
z_disp_wrap,z_plot,normMag,binMask,binv,rSize,monitorSize,scrnIdx, ...
coordType,abs_ang,MaxMask,plotMask)
% Generates and plots a Laguerre Gauss + Zernike mask
%
% Inputs: 
%  X,Y: A grid of the spatial vector: 2D Cartesian coordiantes. They
%  may already have the shiftCart
%  r,phi: polar coordinates (r in cm)
%  phaseValues: discretized phi vector on [-pi,pi].
%                       gl = length(PhaseValues): number of grey levels 
%  tc: Topological charge
%  s: Sign of mask (+1 or -1)
%  ph0: initial phase of the spiral phase mask
%  p: number of radial nodes. If p=0, normal helicoid masks are obtained.
%     If they are used and tc=0(m=0); binary masks are obtained.
%     Even p; rings are ones. Odd p; rings are zeroes. Use mask = mask'
%  WsizeRatio: width of the modes: related with the radius of the phase and
%              with the disks on the magnitude
%  z_coeff: 
%     selected aberrations: any combination on [1,20]; zernike
%     coefficients. ANSI standard. Zernike with desired weights.
%     Put the aberration j-numbers: j = [0 1 2 3 4 ... 20]
%     A) Random surface: -1; Zernike with random weights on [1,20]
%     B) Column vector with normalized coefficients on [0,0.99]
%        [NOll weight]
%     C) Row vector with the j #s that one wants to plot (only integers)
%        [Noll indices]
%
%     Examples:
%     A) z_coeff = -1 (just that)
%     B) z_coeff = [0.1 -0.13 0.15 0 -0.99 0.78] 
%     C) z_coeff = [2 3 1 5 7] (normal)
%  z_a: arbitrary constant; the bigger, the more intense the mask;ref: a=20
%  L: laser wavelength [um]
%  z_pupil: defines pupil relative size (w.r.t. sSize), like a percentage
%  disp_wrap: original (0) or wrapped mask (1)
%  plot_z: plot (1); no plot (0)
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
%  abs_ang: custom(0)[mask real-valued]; magnitude (1); phase (2)
%  MaxMask: defines if the mask should be maximized when coordType = 1
%           -0: custom-size mask that depends on the variable sSize   
%           -1: maximizes the mask for coordType = 1
%           -2: maximized mask but keeping its rectangular fashion
%  plotMask:  no (0); on the screen (1); on the SLM (2); on the screen, but
%             a surface (3)
%
% Outputs:
% mask: LG+Zernike mask. Complex structure that has not been truncated and
% is wrapped on [-pi,pi]. mask = exp(i*UnwrappedMask).
%  wrapMask: truncation and angle operations on mask.
%  wrapMaskFig: figure handler if needed outside the function

showEachMask = 0; % 0: so that not all of the individual masks are shown
                  % during the measurements. Same as plotMask
MaskPupil = 0; % Always 0, since Zernike has its own pupil
[maskZ,~,~]= f_ZernikeMask(X,Y,r,z_coeff,z_a,L,phaseValues, ...
z_pupil,z_disp_wrap,z_plot,normMag,binMask,binv,rSize,monitorSize, ...
scrnIdx,coordType,MaxMask,showEachMask);
[maskLG,~,~] = f_LGMask(r,phi,phaseValues,tc,s,ph0,p,WsizeRatio, ...
normMag,binMask,binv,MaskPupil,rSize,monitorSize,scrnIdx,coordType, ...
abs_ang,MaxMask,showEachMask);

%% LG mask matrix truncation so that its size fits maskZ
A = size(maskLG) - size(maskZ); % Size comparison
% Sizes differ when coordType = 2
if any(A) % True when tests whether any of the elements along various
          % dimensions of an array are nonzero
    [ZmidX,ZmidY] = f_ComputeMatrixMidPoints(maskZ);
    [LGmidX,LGmidY] = f_ComputeMatrixMidPoints(maskLG);
    maskLG = f_TruncateMatrix(maskZ,ZmidX,ZmidY,maskLG,LGmidX,LGmidY);
end 

%% Point-to-point masks product    
mask = maskZ.*maskLG; % Combined mask. It is recommended to use binary
                      % masks on LG
  
%% Plot the combined LG + Z mask 
gl = length(phaseValues); % Number of grey levels 
tit = strcat('LG phase mask with tc=',num2str(tc), ...
             {' '},', p=',num2str(p),{' '}, ...
             '+ Zernike mask',{' '},'with',{' '},'gl=',num2str(gl));
str = ''; % Empty, it only works for abs_ang = 0
[wrapMask,wrapMaskFig] = f_ProjectMask(r,mask,phaseValues,normMag, ...
binMask,binv,MaskPupil,rSize,monitorSize,scrnIdx,tit,str,coordType, ...
abs_ang,MaxMask,plotMask);
end
