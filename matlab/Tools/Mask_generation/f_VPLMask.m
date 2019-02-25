%% Vortex-Producing Lens (VPL) phase mask
% Taken from
% 1_edgar_2013_High-quality optical vortex-beam generation_E-Rueda_OL.pdf
% Equation 3, page 2

function mask = f_VPLMask(r,phi,gl,glphi,mingl,maxgl,levShft,tc,s,ph0, ...
L,f_FR,normMag,binMask,binv,monitorSize,scrnIdx,coordType,abs_ang,plotMask)
% Generates and plots a VPL mask:  helicoidal mask + fresnel lens
%
% Inputs: 
%  x,y: cartesian coordinates vector
%  r,phi: polar coordinates (r in cm)
%  gl: number of grey levels (normally 256)
%  glphi: discretized phi vector on [-pi,pi].
%  mingl,maxgl: minimum/maximum gray level depth. Ref: 0,255
%  levShft: corresponds to the brightness or constant shift of the gl's
%  tc: Topological charge
%  s: Sign of mask (+1 or -1)
%  ph0: initial phase of the spiral phase mask
%  L: Laser wavelength [um]
%  f_FR: Fresnel lens focal distance or diffractive lens phase focal length
%        in um
%  binMask: binarizes the mask w.r.t the max and min of the phase (boolean)
%  showM: show the mask. yes(1); no(0)
%
% Outputs:
% mask: Vortex-Producing Lens (VPL) phase mask. Complex structure that has
% not been truncated and is wrapped on [-pi,pi]. mask = exp(i*UnwrppedMsk).
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
f_ProjectMask(r,mask,gl,glphi,mingl,maxgl,levShft,normMag,binMask, ...
binv,monitorSize,scrnIdx,tit,str,coordType,abs_ang,plotMask);

end