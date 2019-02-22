%% Spiral Phase Mask
function mask = f_Spiral_Mask(x,y,r,phi,gl,glphi,mingl,maxgl,levShft,tc,s,ph0,binMask,showM)
% Plots a custom spiral phase mask with a specific topological charge
% and an initial angle. Can be plotted on the SLM screen or normally
%
% Inputs: 
%  x,y: cartesian coordinates vector
%  phi,r: polar coordinates (r in cm)
%  gl: number of grey levels (normally 256)
%  glphi: discretized phi vector on [-pi,pi].
%  mingl,maxgl: minimum/maximum gray level depth. Ref: 0,255
%  levShft: corresponds to the brightness or constant shift of the gl's
%  tc: Topological charge
%  s: Sign of mask (+1 or -1)
%  ph0: initial phase of the spiral phase mask
%  binMask: binarizes the mask w.r.t the max and min of the phase (boolean)
%  showM: show the mask. yes(1); no(0)
%
% Outputs:
%  mask: spiral phase mask

%% Spiral phase mask Generation
% m = -s*tc; % OLD: Minus to compensate convention
m = s*tc; % tc with a sign
mask = m*(phi + ph0); % General mask. Angle phi is wrapped on [-pi,pi]

%% Discretized phase mask
% glPHI = meshgrid(glphi);
% discretPhi = f_discretizeMask(mask,glphi); % Mask discretization
% mask = f_scaleMatrix(discretPhi,mingl,maxgl) + levShft; 

mask = exp(1i*mask); % Wrapped mask and complex
% mask = angle(mask); % OLD

% OLD:
%  mask = angle(exp(1i*m*(phi+ph0+pi))) + pi; % Wrapped on [0,2*pi]
%  mask = mod(m*(phi+ph0),2*pi);
%  pi was added to compensate the initial angle and to correct the wrapping
%  interval

%% Circular pupil and wrapping
wrappedMask = f_mask_circ_angle_gl(r,mask,binMask,glphi,mingl,maxgl,levShft);

%% Plot (with axes)
tit = strcat('Spiral phase mask with topological charge ',num2str(tc));
f_fig_maskPCscreen(x, y, wrappedMask, tit, gl, showM);

end