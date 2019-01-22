%% Laguerre Gauss phase masks

function mask = f_LG_Mask(x,y,r,phi,gl,tc,s,ph0,p,W,binv,norm,abs_ang,binMask,showM)
% Inputs: 
%  x,y: cartesian coordinates
%  r,phi: polar coordinates (r in cm)
%  gl: number of grey levels (normally 256)
%  tc: topological charge = m (azimuthal index)
%  s: sign of mask (+1 or -1)
%  ph0: initial phase of the azimuthal part (and of the spiral phase mask)
%  p: number of radial nodes. If p=0, normal helicoid masks are obtained.
%     If they are used and tc=0(m=0); binary masks are obtained.
%     Even p; rings are ones. Odd p; rings are zeroes. Use mask = mask'
%  W: width of the modes: related with the radius of the phase and with the
%     disks on the magnitude
%  bininv: Binary inversion of the mask. Only applied when tc=0. Binary 
%          masks are only abtained when tc=0. yes(1); no(0)
%  norm: normalize magnitude and phase. yes(1); no(0)
%  abs_ang: Magnitude (1); Phase (2)
%  binMask: binarizes the mask w.r.t the max and min of the phase (boolean)
%  showM: show the mask. yes(1); no(0)
%
% Output:
%  mask: Phase mask
%
% Samuel Plazas Escudero - Juan José Cadavid - 2018 - Advanced Project 1

%% Parameters: Laguerre-Gauss
m = tc; % Azimuthal index = topological charge
sSize = length(x); % Size of the original x,y coordinates
W = W/sSize; % Normalization with number of samples

%% Phase mask
LG = f_LaguerreGauss(r,phi,m,s,ph0,p,W); % Generate Laguerre-Gauss mode: it
                                         % has both magnitude and phase
mag = abs(LG); % Magnitude of LG
% mask = mag.*exp(1i*angle(LG)); % Redundant step
mask = LG; % Phase mask. Wrapped on [-pi,pi], modulo(2pi)

%% Binary inversion
if binv == 1 && tc == 0 % Only applies binary inversion when tc = 0
    mask = ctranspose(-mask); % Binary inversion 
end

%% Normalization constants (amplitude and phase)
wrappedMask = f_circularPupil_maskAngle(r,mask,binMask);

if norm == 1
    norm_ang = max(max(wrappedMask)); % Max value
    wrappedMask = wrappedMask/norm_ang; % Normalization
    norm_mag = max(max(mag)); % Max value
    mag = mag/norm_mag; % Normalization
end

%% Plot binary mask
if showM == 1
  if abs_ang == 2
    tit = ['LG phase mask with topological charge ' num2str(tc) ...
           ' and radial node ' num2str(p)];   
    f_fig_maskPCscreen(x, y, wrappedMask, tit, gl, showM);
  else % abs_ang == 1
    plotMask = showM; % plotMask = show; for 0 and 1.
    f_fig_maskSLM(x,y,r,mag,0,0,0,0,gl,abs_ang,binMask,plotMask); 
    title('Amplitude of LG');
    cbh = colorbar; cbh.Label.String = 'Value';
  end
end
% f_fig_maskSLM(x,y,r,mask,m,n,a,b,gl,abs_ang,binMask,plotMask)

%%% OLD
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
