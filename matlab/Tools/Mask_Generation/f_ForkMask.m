%% Elliptic Gaussian Vortex (EGV)pi
% Taken from: 2015_Vortex_CGH_Adjustable-SPP_Jain

function [mask,wrapMask,wrapMaskFig] = f_ForkMask(X,Y,r,phi,phaseValues,...
tc,s,ph0,period,T0,frkTyp,Aalpha,Angalp,Angbet,normMag,binMask,binv, ...
MaskPupil,rSize,monitorSize,scrnIdx,coordType,abs_ang,MaxMask,plotMask)
% Plots a custom spiral phase mask with a specific topological charge
% and an initial angle. Can be plotted on the SLM screen or normally
%
% Inputs: 
%  X,Y: A grid of the spatial vector: 2D Cartesian coordiantes
%  r,phi: polar coordinates (r in cm)
%  phaseValues: discretized phi vector on [-pi,pi].
%                       gl = length(PhaseValues): number of grey levels 
%  tc: Topological charge
%  s: Sign of mask (+1 or -1)
%  ph0: initial phase of the spiral phase mask
%  period: of the grating (fringe spacing)
%  T0: constant absorption coefficient of the hologram
%  frkTyp: 1: smooth transition; 2: phase jump transition
%  Aalpha: amplitude of the phase modulation
%  Angalp,Angbet: diffraction angles of horizontal and vertical directions
%                 they are limited to [-pi/2,pi/2] [radians]
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
%  abs_ang: custom(0)[str has to be defined for this case], magnitude
%           (1) or phase (2) plot. Doesn't apply for Zernike and LG +
%           Zernike.
%  MaxMask: defines if the mask should be maximized when coordType = 1
%           -0: custom-size mask that depends on the variable sSize   
%           -1: maximizes the mask for coordType = 1
%           -2: maximized mask but keeping its rectangular fashion
%  plotMask:  no (0); on the screen (1); on the SLM (2); on the screen, but
%             a surface (3)
%
% Outputs: Threefold dislocation hologram or double pitch fork hologram or
%          fork phase mask or holograma en forma de tenedor
% mask: Fork mask. Complex structure that has not been truncated and is
%       wrapped on [-pi,pi]. mask = exp(i*UnwrappedMask).
%  wrapMask: truncation and angle operations on mask.
%  wrapMaskFig: figure handler if needed outside the function

switch frkTyp
    case 1 % Smooth transition hologram
        %% Double pitch fork hologram phase
        maskFork = -(2*pi/period)*r.*cos(phi);

        %% Spiral phase mask Generation
        m = s*tc; % tc with a sign
        maskSPP = m*(phi + ph0); % General mask. Angle phi is wrapped on [-pi,pi]

        %% Fork Mask
        mask = T0*exp(1i*Aalpha*cos(maskSPP + maskFork)); % Wrapped mask

    case 2 % Phase jump hologram
        %% Spiral phase mask Generation
        m = s*tc; % tc with a sign
        maskSPP = m*(phi + ph0); % General mask. Angle phi is wrapped on [-pi,pi]
        
        %% Double pitch fork hologram phase
        maskFork = (2*pi/period)*(sin(Angalp*X) + sin(Angbet*Y));
        % Original: (2*pi/L) as it should approximate the wavelength to
        % have considerable effects of diffraction
        % Test of the unwrapped fork:
        % figure, imagesc(maskFork);
        
        %% Fork Mask
        mask = exp(1i*(maskFork+maskSPP));      
       
end

%% Plot the mask
gl = length(phaseValues); % Number of grey levels 
tit = strcat('Fork mask with tc=',num2str(tc), ...
             ', period=',num2str(period),{' '},'and',{' '},'gl=',num2str(gl));
str = ''; % Empty, it only works for abs_ang = 0
[wrapMask,wrapMaskFig] = f_ProjectMask(r,mask,phaseValues,normMag, ...
binMask,binv,MaskPupil,rSize,monitorSize,scrnIdx,tit,str,coordType, ...
abs_ang,MaxMask,plotMask);

end