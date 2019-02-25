%% Elliptic Gaussian Vortex (EGV)pi
% Taken from: 2015_Vortex_CGH_Adjustable-SPP_Jain

function mask = f_ForkMask(X,Y,r,phi,gl,glphi,mingl,maxgl,levShft,tc,s, ...
ph0,L,period,T0,frkTyp,Aalpha,Angalp,Angbet,normMag,binMask,binv, ...
monitorSize,scrnIdx,coordType,abs_ang,plotMask)
% Plots a custom spiral phase mask with a specific topological charge
% and an initial angle. Can be plotted on the SLM screen or normally
%
% Inputs: 
%  X,Y: A grid of the spatial vector: 2D Cartesian coordiantes
%  r,phi: polar coordinates (r in cm)
%  gl: number of grey levels (normally 256)
%  glphi: discretized phi vector on [-pi,pi].
%  mingl,maxgl: minimum/maximum gray level depth. Ref: 0,255
%  levShft: corresponds to the brightness or constant shift of the gl's
%  tc: Topological charge
%  s: Sign of mask (+1 or -1)
%  ph0: initial phase of the spiral phase mask
%  L: Laser wavelength [um]
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
%  monitorSize: size of the selected screen 
%  scrnIdx: screen number selector. In [1,N] with N the # of screen
%  coordType: type of calculation of the spatial coordinates. def: 2 
%    -1: size defined by the user, space support defined by the SLM to use
%    -2: size defined by the resolution of the selected screen 
%  abs_ang: custom(0)[str has to be defined for this case], magnitude
%           (1) or phase (2) plot. Doesn't apply for Zernike and LG +
%           Zernike.
%  plotMask:  no (0); on the screen (1); on the SLM (2); on the screen, but
%             a surface (3)
%
% Outputs: Threefold dislocation hologram or double pitch fork hologram or
%          fork phase mask or holograma en forma de tenedor
% mask: Fork mask. Complex structure that has not been truncated and is
%       wrapped on [-pi,pi]. mask = exp(i*UnwrappedMask).

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
        % Original: (2*pi/L)
        
        %% Fork Mask
        mask = exp(1i*(maskFork+maskSPP));
       
        % Test of the unwrapped fork
        % figure, imagesc(maskFork);
        
end

%% Plot the mask
tit = strcat('Fork mask with topological charge',{' '},num2str(tc), ...
             {' '},'and period =',{' '},num2str(period));
str = ''; % Empty, it only works for abs_ang = 0
f_ProjectMask(r,mask,gl,glphi,mingl,maxgl,levShft,normMag,binMask, ...
binv,monitorSize,scrnIdx,tit,str,coordType,abs_ang,plotMask);

end