%% Laguerre Gauss phase masks

function mask = f_LGMask(r,phi,gl,glphi,mingl,maxgl,levShft,tc,s,ph0,p, ...
W,normMag,binMask,binv,monitorSize,scrnIdx,coordType,abs_ang,plotMask)
% Inputs: 
%  r,phi: polar coordinates (r in cm) for both the PC and SLM
%  gl: number of grey levels (normally 256)
%  glphi: discretized phi vector on [-pi,pi].
%  mingl,maxgl: minimum/maximum gray level depth. Ref: 0,255
%  levShft: corresponds to the brightness or constant shift of the gl's
%  tc: topological charge = m (azimuthal index). When tc = 0, the mask is
%      purely real
%  s: sign of mask (+1 or -1)
%  ph0: initial phase of the azimuthal part (and of the spiral phase mask)
%  p: number of radial nodes. If p=0, normal helicoid masks are obtained.
%     If they are used and tc=0(m=0); binary masks are obtained.
%     Even p; rings are ones. Odd p; rings are zeroes. Use mask = mask'
%  W: width of the modes: related with the radius of the phase and with the
%     disks on the magnitude
%  normMag: normalize magnitude. yes(1); no(0)
%  binMask: binarizes the mask w.r.t the max and min of the phase (boolean)
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
%  mask: LG mask. Complex structure that has not been truncated and is 
%  wrapped on [-pi,pi]. mask = exp(i*UnwrappedMask).
%
% Samuel Plazas Escudero - Juan Jos� Cadavid - 2018 - Advanced Project 1

%% Parameters: Laguerre-Gauss
m = tc; % Azimuthal index = topological charge
% sSize = length(x); % Size of the original x,y coordinates
% sSize = min(size(r)); % Addded but not fully sure
W = W/100; % Normalization with number of samples

%% Laguerre-Gauss mask
mask = f_LaguerreGauss(r,phi,m,s,ph0,p,W); % Generates a Laguerre-Gauss 
% mode: it has both magnitude and phase, meaning that:
% mask = abs(mask).*exp(1i*angle(mask))

% OLD:
% mag = abs(LG); % Magnitude of LG
% mask = LG; % Phase mask. Wrapped on [-pi,pi], modulo(2pi)

%% Plot the mask
tit = strcat('LG phase mask with topological charge',{' '},num2str(tc), ...
                 ' and radial node',{' '},num2str(p));   
f_ProjectMask(r,mask,gl,glphi,mingl,maxgl,levShft,normMag,binMask,binv, ...
              monitorSize,scrnIdx,tit,coordType,abs_ang,plotMask);
          
%%% OLD
%  title('Amplitude of LG');
%  cbh = colorbar; cbh.Label.String = 'Value';
% if abs_ang == 2
%   h = pcolor(x,y,mask); 
% else % var == 1
%   h = pcolor(x,y,abs(LG));
%   % Only positive values since it is a intensity measurement
% end
% colormap('bone'), set(h,'EdgeColor','none'), set(h,'FaceColor','interp');
% set(gca,'Visible','off'), set(gcf,'Color','black'), axis square, hold off;
% shg; % shg makes the current figure visible and raises it above all other
%      % figures on the screen.

end
