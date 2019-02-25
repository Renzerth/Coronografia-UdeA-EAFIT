%% Spiral Phase Mask
function [mask,wrapMask,wrapMaskFig] = f_SpiralMask(r,phi,gl,phaseValues, ...
mingl,maxgl,levShft,tc,s,ph0,normMag,binMask,binv,monitorSize,scrnIdx, ...
coordType,abs_ang,plotMask)
% Plots a custom spiral phase mask with a specific topological charge
% and an initial angle. Can be plotted on the SLM screen or normally
%
% Inputs:
%  r,phi: polar coordinates for both the PC and SLM
%  gl: number of grey levels (normally 256)
%  phaseValues: discretized phi vector on [-pi,pi].
%  mingl,maxgl: minimum/maximum gray level depth. Ref: 0,255
%  levShft: corresponds to the brightness or constant shift of the gl's
%  tc: Topological charge
%  s: Sign of mask (+1 or -1)
%  ph0: initial phase of the spiral phase mask
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
%  plotMask:  no (0); on the screen (1); on the SLM (2); on the screen, but
%             a surface (3)
%
% Outputs:
%  mask: spiral phase mask. Complex structure that has not been truncated 
%  and is wrapped on [-pi,pi]. mask = exp(i*UnwrappedMask).
%  wrapMask: truncation and angle operations on mask.
%  wrapMaskFig: figure handler if needed outside the function

%% Spiral phase mask generation
m = s*tc; % tc with a sign
mask = m*(phi + ph0); % General mask. Angle phi is wrapped on [-pi,pi]
mask = exp(1i*mask); % Wrapped mask and complex
% mask = mod(mask,2*pi) - pi; % Equivalent operation

%% Plot the mask
tit = strcat('Spiral phase mask with topological charge',{' '}, ...
      num2str(tc));
str = ''; % Empty, it only works for abs_ang = 0
[wrapMask,wrapMaskFig] = f_ProjectMask(r,mask,gl,phaseValues,mingl,maxgl, ...
levShft,normMag,binMask,binv,monitorSize,scrnIdx,tit,str,coordType, ...
abs_ang,plotMask);

end