%% Elliptic Gaussian Vortex (EGV)
% Taken from: 1_2017_IMP_Elliptic_Gaussian Optical Vortices_PRA_Kotlyar

function mask = f_EGV_Mask(x,y,X,Y,r,gl,tc,s,ph0,bcst,binMask,showM)
% Plots a custom spiral phase mask with a specific topological charge
% and an initial angle. Can be plotted on the SLM screen or normally
%
% Inputs: 
%  x,y: cartesian coordinates vector
%  X,Y: A symmetric grid of the spatial vector: 2D Cartesian coordiantes
%  r: polar coordinate (in cm)
%  gl: number of grey levels (normally 256)
%  tc: Topological charge
%  s: Sign of mask (+1 or -1)
%  ph0: initial phase of the spiral phase mask
%  bcst (beta): Ellipticity.cy/cx = 1/alpha. Ref: .1, .2, .4, .6, .8 and 1
%  alpha is a dimensionless parameter that defines the ellipticity of
%  the intensity null: if alpha < 1, the major axis is on the x axis,
%  if alpha > 1 – on the y axis, if alpha < 0, the vortex phase rotates
%  clockwise, if alpha > 0 – anticlockwise.
%  binMask: binarizes the mask w.r.t the max and min of the phase (boolean)
%  showM: show the mask. yes(1); no(0)
%
% Outputs:
% mask: Elliptic Gaussian Vortex (EGV)
%


%% Elliptic Gaussian beam phase mask
phi = atan2(bcst*Y,X); % Angle component with elliptic gaussian beam
% atan2(beta^(0.5)*Y,beta^-(0.5)*X) % OLD: Works the same as above

%% Spiral phase mask Generation
m = s*tc; % tc with a sign
mask = m*(phi + ph0); % General mask. Angle phi is wrapped on [-pi,pi]
mask = exp(1i*mask); % Wrapped mask

%% Plot (with axes)
if showM == 1
  wrappedMask = f_circularPupil_maskAngle(r,mask,binMask);
  tit = ['EGV with topological charge ' num2str(tc) ...
         ' and beta = ' num2str(bcst)];
  f_fig_maskPCscreen(x, y, wrappedMask, tit, gl, showM);
end

end