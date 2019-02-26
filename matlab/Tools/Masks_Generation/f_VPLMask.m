%% Vortex-Producing Lens (VPL) phase mask
% Taken from
% 1_edgar_2013_High-quality optical vortex-beam generation_E-Rueda_OL.pdf
% Equation 3, page 2

function [mask,wrapMask,wrapMaskFig] = f_VPLMask(r,phi,gl,phaseValues,mingl, ...
maxgl,levShft,tc,s,ph0,L,f_FR,normMag,binMask,binv,monitorSize,scrnIdx, ...
coordType,abs_ang,MaxMask,plotMask)
% Generates and plots a VPL mask:  helicoidal mask + fresnel lens
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
%  L: Laser wavelength [um]
%  f_FR: Fresnel lens focal distance or diffractive lens phase focal length
%        in um
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
%  MaxMask: maximizes the mask for coordType = 1 (0): doesn't
%  plotMask:  no (0); on the screen (1); on the SLM (2); on the screen, but
%             a surface (3)
%
% Outputs:
%  mask: Vortex-Producing Lens (VPL) phase mask. Complex structure that has 
%  not been truncated and is wrapped on [-pi,pi]. mask=exp(i*UnwrappedMask)
%  wrapMask: truncation and angle operations on mask.
%  wrapMaskFig: figure handler if needed outside the function
%

%% Spiral phase mask Generation
m = s*tc;  % tc with a sign
maskSPP = m*(phi + ph0); % General mask. Angle phi is wrapped on [-pi,pi]
maskSPP = exp(1i*maskSPP); % Wrapped mask

%% VPL Mask generation
scaleFactor = 1e4; % cm to um. Constant factor
r = r*scaleFactor; % r converted to um
VPLfactor = -(pi*r.^2)/(L*f_FR);
maskVPL = exp(1i*VPLfactor);
mask = maskSPP.*maskVPL; 
% maskSPP*maskVPL; gives cool results as the Fresnel lens but in a linear
% fashion: a chirp ramp for tc = 0. It is just a test without validity

%% Plot the mask
tit = strcat('VPL with topological charge',{' '},num2str(tc),{' '}, ...
             'and',{' '},num2str(gl),{' '},'gray levels');  
str = ''; % Empty, it only works for abs_ang = 0
[wrapMask,wrapMaskFig] = f_ProjectMask(r,mask,gl,phaseValues,mingl,maxgl, ...
levShft,normMag,binMask,binv,monitorSize,scrnIdx,tit,str,coordType, ...
abs_ang,MaxMask,plotMask);

end