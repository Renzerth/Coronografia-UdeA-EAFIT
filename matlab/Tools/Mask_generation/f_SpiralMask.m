%% Spiral Phase Mask
function mask = f_SpiralMask(r,phi,gl,glphi,mingl,maxgl,levShft,tc,s, ...
ph0,normMag,binMask,binv,monitorSize,scrnIdx,coordType,abs_ang,plotMask)
% Plots a custom spiral phase mask with a specific topological charge
% and an initial angle. Can be plotted on the SLM screen or normally
%
% Inputs:
%  r,phi: polar coordinates (r in cm) for both the PC and SLM
%  gl: number of grey levels (normally 256)
%  glphi: discretized phi vector on [-pi,pi].
%  mingl,maxgl: minimum/maximum gray level depth. Ref: 0,255
%  levShft: corresponds to the brightness or constant shift of the gl's
%  tc: Topological charge
%  s: Sign of mask (+1 or -1)
%  ph0: initial phase of the spiral phase mask
%  normMag: normalize magnitude. yes(1); no(0)
%  binMask: binarizes the mask w.r.t the max and min of the phase (boolean)
%  monitorSize: size of the selected screen 
%  binv: binary inversion of the mask: yes(1); no(0). Only applies when 
%        binMask=1. It is usefull to be applied for odd p's on LG beams
%  monitorSize: size of the selected screen 
%  binv: binary inversion of the mask: yes(1); no(0). Only applies when 
%        binMask=1. It is usefull to be applied for odd p's on LG beams
%  scrnIdx: screen number selector. In [1,N] with N the # of screen
%  coordType: type of calculation of the spatial coordinates. def: 2 
%    -1: size defined by the user, space support defined by the SLM to use
%    -2: size defined by the resolution of the selected screen    
%  abs_ang: custom(0)[mask real-valued]; magnitude (1); phase (2)
%  plotMask:  no (0); on the screen (1); on the SLM (2); on the screen, but
%             a surface (3)
%
% Output:
%  mask: spiral phase mask. Complex structure that has not been truncated
%        and is wrapped on [-pi,pi]. mask = exp(i*UnwrappedMask).

%% Spiral phase mask generation
m = s*tc; % tc with a sign
mask = m*(phi + ph0); % General mask. Angle phi is wrapped on [-pi,pi]
mask = exp(1i*mask); % Wrapped mask and complex

% OLD:
% m = -s*tc; %  Minus to compensate convention

% OLD1:
% glPHI = meshgrid(glphi);
% discretPhi = f_discretizeMask(mask,glphi); % Mask discretization
% mask = f_scaleMatrix(discretPhi,mingl,maxgl) + levShft; 

% OLD2:
% mask = angle(mask);
%  mask = angle(exp(1i*m*(phi+ph0+pi))) + pi; % Wrapped on [0,2*pi]
%  mask = mod(m*(phi+ph0),2*pi);
%  pi was added to compensate the initial angle and to correct the wrapping
%  interval

%% Plot the mask
tit = strcat('Spiral phase mask with topological charge',{' '}, ...
      num2str(tc));
f_ProjectMask(r,mask,gl,glphi,mingl,maxgl,levShft,normMag,binMask,binv, ...
              monitorSize,scrnIdx,tit,coordType,abs_ang,plotMask);

end