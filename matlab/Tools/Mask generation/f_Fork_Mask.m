%% Elliptic Gaussian Vortex (EGV)pi
% Taken from: 2015_Vortex_CGH_Adjustable-SPP_Jain

function mask = f_Fork_Mask(x,y,X,Y,r,phi,gl,tc,s,ph0,L,period,T0,frkTyp,Aalpha,Angalp,Angbet,binMask,showM)
% Plots a custom spiral phase mask with a specific topological charge
% and an initial angle. Can be plotted on the SLM screen or normally
%
% Inputs: 
%  x,y: cartesian coordinates vector
%  X,Y: A symmetric grid of the spatial vector: 2D Cartesian coordiantes
%  phi,r: polar coordinates (r in cm)
%  gl: number of grey levels (normally 256)
%  tc: Topological charge
%  s: Sign of mask (+1 or -1)
%  ph0: initial phase of the spiral phase mask
%  L: Laser wavelength [um]
%  period: of the grating (fringe spacing)
%  T0: constant absorption coefficient of the hologram
%  frkTyp: 1: smooth transition; 2: phase jump transition
%  Aalpha: amplitude of the phase modulation
%  binMask: binarizes the mask w.r.t the max and min of the phase (boolean)
%  Angalp,Angbet: diffraction angles of horizontal and vertical directions
%                 they are limited to [-pi/2,pi/2] [radians]
%  showM: show the mask. yes(1); no(0)
%
% Outputs: Threefold dislocation hologram or double pitch fork hologram or
%          fork phase mask or holograma en forma de tenedor
% mask: Elliptic Gaussian Vortex (EGV)
%

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
        
end

%% Plot (with axes)
if showM == 1
  wrappedMask = f_circularPupil_maskAngle(r,mask,binMask);
  tit = ['Fork mask with topological charge ' num2str(tc) ...
         ' and period = ' num2str(period)];
  f_fig_maskPCscreen(x, y, wrappedMask, tit, gl, showM);
end

end