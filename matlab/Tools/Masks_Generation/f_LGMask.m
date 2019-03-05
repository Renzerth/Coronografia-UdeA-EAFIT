%% Laguerre Gauss phase masks

function [mask,wrapMask,wrapMaskFig] = f_LGMask(r,phi,gl,phaseValues,mingl, ...
maxgl,levShft,tc,s,ph0,p,WsizeRatio,normMag,binMask,binv,MaskPupil,rSize,monitorSize,scrnIdx, ...
coordType,abs_ang,MaxMask,plotMask)
% Inputs: 
%  r,phi: polar coordinates for both the PC and SLM
%  gl: number of grey levels (normally 256)
%  phaseValues: discretized phi vector on [-pi,pi].
%  mingl,maxgl: minimum/maximum gray level depth. Ref: 0,255
%  levShft: corresponds to the brightness or constant shift of the gl's
%  tc: topological charge = m (azimuthal index). When tc = 0, the mask is
%      purely real
%  s: sign of mask (+1 or -1)
%  ph0: initial phase of the azimuthal part (and of the spiral phase mask)
%  p: number of radial nodes. If p=0, normal helicoid masks are obtained.
%     If they are used and tc=0(m=0); binary masks are obtained.
%     Even p; rings are ones. Odd p; rings are zeroes. Use mask = mask'
%  WsizeRatio: width of the modes: related with the radius of the phase and with the
%     disks on the magnitude
%  normMag: normalize magnitude. yes(1); no(0)
%  binMask: binarizes the mask w.r.t the max and min of the phase (boolean)
%  binv: binary inversion of the mask: yes(1); no(0). Only applies when 
%        binMask=1. It is usefull to be applied for odd p's on LG beams
%  MaskPupil: applies a pupil truncation to the mask: (0): no; (1): yes
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
%  mask: LG mask. Complex structure that has not been truncated and is 
%  wrapped on [-pi,pi]. mask = exp(i*UnwrappedMask).
%  wrapMask: truncation and angle operations on mask.
%  wrapMaskFig: figure handler if needed outside the function
%
% Samuel Plazas Escudero - Juan José Cadavid - 2018 - Advanced Project 1

%% Parameters: Laguerre-Gauss
m = tc; % Azimuthal index = topological charge
sSize = min(size(r)); % Addded but not fully sure
WsizeRatio = sSize*WsizeRatio; % Normalization with number of samples

%% Laguerre-Gauss mask
mask = f_LaguerreGauss(r,phi,m,s,ph0,p,WsizeRatio); % Generates a Laguerre-Gauss 
% mode: it has both magnitude and phase, meaning that:
% mask = abs(mask).*exp(1i*angle(mask))

%% Plot the mask
tit = strcat('LG phase mask with topological charge',{' '},num2str(tc), ...
                 ' and radial node',{' '},num2str(p));   
str = ''; % Empty, it only works for abs_ang = 0
[wrapMask,wrapMaskFig] = f_ProjectMask(r,mask,gl,phaseValues,mingl,...
maxgl,levShft,normMag,binMask,binv,MaskPupil,rSize,monitorSize,scrnIdx,tit, ...
str,coordType,abs_ang,MaxMask,plotMask);

%% Mask in a bone colormap
%   h = pcolor(x,y,mask); 
% colormap('bone'), set(h,'EdgeColor','none'), set(h,'FaceColor','interp');
% set(gca,'Visible','off'), set(gcf,'Color','black'), axis square, hold off;
% shg; % shg makes the current figure visible and raises it above all other
%      % figures on the screen.

end
