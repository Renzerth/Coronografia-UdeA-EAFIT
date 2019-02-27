%% Elliptic Gaussian Vortex (EGV)
% Taken from: 1_2017_IMP_Elliptic_Gaussian Optical Vortices_PRA_Kotlyar

function [mask,wrapMask,wrapMaskFig] = f_EGVMask(X,Y,r,gl,phaseValues,mingl, ...
maxgl,levShft,tc,s,ph0,bcst,normMag,binMask,binv,monitorSize,scrnIdx, ...
coordType,abs_ang,MaxMask,plotMask)
% Plots a custom spiral phase mask with a specific topological charge
% and an initial angle. Can be plotted on the SLM screen or normally
%
% Inputs: 
%  X,Y: A meshgrid of the spatial vector: 2D Cartesian coordinates
%  r: polar coordinate (in cm)
%  gl: number of grey levels (normally 256)
%  phaseValues: discretized phi vector on [-pi,pi].
%  mingl,maxgl: minimum/maximum gray level depth. Ref: 0,255
%  levShft: corresponds to the brightness or constant shift of the gl's
%  tc: Topological charge
%  s: Sign of mask (+1 or -1)
%  ph0: initial phase of the spiral phase mask
%  bcst (beta): Ellipticity.cy/cx = 1/alpha. Ref: .1, .2, .4, .6, .8 and 1
%               alpha is a dimensionless parameter that defines the 
%               ellipticity of the intensity null: if alpha < 1, the major 
%               axis is on the x axis, if alpha > 1 – on the y axis, if 
%               alpha < 0, the vortex phase rotates clockwise, if alpha > 0
%               then anticlockwise.
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
%  mask: Elliptic Gaussian Vortex (EGV). Complex structure that has not 
%  been truncated and is wrapped on [-pi,pi]. mask = exp(i*UnwrappedMask).
%  wrapMask: truncation and angle operations on mask.
%  wrapMaskFig: figure handler if needed outside the function

%% Elliptic Gaussian beam phase mask: scalated azimuthal coordinate
phi = atan2(bcst*Y,X); % Angle component with elliptic gaussian beam
% atan2(beta^(0.5)*Y,beta^-(0.5)*X) % Works the same as above

%% Spiral phase mask Generation
m = s*tc; % tc with a sign
mask = m*(phi + ph0); % General mask. Angle phi is wrapped on [-pi,pi]
mask = exp(1i*mask); % Wrapped mask

%% Plot the mask
tit = strcat('EGV with topological charge',{' '},num2str(tc),{' '}, ...
             'and beta = ',{' '},num2str(bcst));
str = ''; % Empty, it only works for abs_ang = 0
[wrapMask,wrapMaskFig] = f_ProjectMask(r,mask,gl,phaseValues,mingl,maxgl, ...
levShft,normMag,binMask,binv,monitorSize,scrnIdx,tit,str,coordType, ...
abs_ang,MaxMask,plotMask);

end