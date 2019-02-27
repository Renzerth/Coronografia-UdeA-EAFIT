%% Laguerre Gauss Binary masks + Zernike Phases
function [mask,wrapMask,wrapMaskFig] = f_LGZernikeMask(r,phi,gl,phaseValues, ...
mingl,maxgl,levShft,tc,s,ph0,p,W,z_coeff,a,frac,L,pupil,ZernikeSize,disp_wrap,...
plot_z,normMag,binMask,binv,monitorSize,scrnIdx,coordType,abs_ang,MaxMask,plotMask)
% Generates and plots a Laguerre Gauss + Zernike mask
%
% Inputs: 
%  r,phi: polar coordinates (r in cm)
%  gl: number of grey levels (normally 256)
%  phaseValues: discretized phi vector on [-pi,pi].
%  mingl,maxgl: minimum/maximum gray level depth. Ref: 0,255
%  levShft: corresponds to the brightness or constant shift of the gl's
%  tc: Topological charge
%  s: Sign of mask (+1 or -1)
%  ph0: initial phase of the spiral phase mask
%  p: number of radial nodes. If p=0, normal helicoid masks are obtained.
%     If they are used and tc=0(m=0); binary masks are obtained.
%     Even p; rings are ones. Odd p; rings are zeroes. Use mask = mask'
%  W: width of the modes: related with the radius of the phase and with the
%     disks on the magnitude
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
%  a: arbitrary constant; the bigger, the more intense; ref: a=20
%  frac: to adjust the wrapped phase; ref: 0.125
%  L: laser wavelength [um]
%  pupil: defines pupil relative size (w.r.t. sSize), like a percentage
%  ZernikeSize: screen size for the Zernike polynomials generation
%  disp_wrap: original (0) or wrapped mask (1)
%  plot_z: plot (1); no plot (0)
%  normMag: normalize magnitude. yes(1); no(0)
%  binMask: binarizes the mask w.r.t the max and min of the phase (boolean)
%  binv: binary inversion of the mask: yes(1); no(0). Only applies when 
%        binMask=1. It is usefull to be applied for odd p's on LG beams
%  monitorSize: size of the selected screen 
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

showEachMask = 0;
[maskZ,~,~]= f_ZernikeMask(r,z_coeff,a,frac,L,gl,phaseValues,mingl,maxgl, ...
levShft,pupil,ZernikeSize,disp_wrap,plot_z,normMag,binMask,binv,monitorSize, ...
scrnIdx,coordType,showEachMask);
[maskLG,~,~] = f_LGMask(r,phi,gl,phaseValues,mingl,maxgl,levShft,tc,s,ph0,p, ...
W,normMag,binMask,binv,monitorSize,scrnIdx,coordType,abs_ang,showEachMask);

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
% maskZ = maskZ*100; % Just tests
mask = maskZ.*maskLG; % Combined mask. It is recommended to use binary
                      % masks on LG
  
%% Plot the combined LG + Z mask                      
tit = strcat('LG phase mask with topological charge',{' '},num2str(tc), ...
             {' '},'and radial node',{' '},num2str(p),{' '}, ...
             '+ Zernike mask');
str = ''; % Empty, it only works for abs_ang = 0         
[wrapMask,wrapMaskFig] = f_ProjectMask(r,mask,gl,phaseValues,mingl,maxgl, ...
levShft,normMag,binMask,binv,monitorSize,scrnIdx,tit,str,coordType, ...
abs_ang,MaxMask,plotMask);
end
