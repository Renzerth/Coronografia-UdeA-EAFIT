%% Generate a Zernike phase mask (wavefront)
function mask = f_Zernike_Mask(x,y,r,z_coeff,a,frac,L,gl,glphi,mingl, ...
                  maxgl,levShft,pupil,sSize,disp_wrap,plot_z,binMask,monitorSize,showM)
% Characterizes the aberrations of the system
% Inputs:
%  x,y: cartesian coordinates vector
%  r: polar coordinate (in cm)
%  z_coeff: 
%     selected aberrations: any combination on [1,20]; zernike
%     coefficients. ANSI standard. Zernike with desired weights.
%     Put the aberration j-numbers: j = [0 1 2 3 4 ... 20]
%     A) Random surface: -1; Zernike with random weights on [1,20]
%     B) Column vector with normalized coefficients on [0,1]
%     C) Row vector with the j #s that one wants to plot (only integers)
%
%     Examples:
%     A) z_coeff = -1 (just that)
%     B) z_coeff = [0.1 -0.13 0.15 0 -1 0.78]' (transposed)
%     C) z_coeff = [2 3 1 5 7] (normal)
%  a: arbitrary constant; the bigger, the more intense; ref: a=20
%  frac: to adjust the wrapped phase; ref: 0.125
%  L: laser wavelength [um]
%  pupil
%  gl: gray levels of Zernike. Normally 256
%  glphi: discretized phi vector on [-pi,pi].
%  mingl,maxgl: minimum/maximum gray level depth. Ref: 0,255
%  levShft: corresponds to the brightness or constant shift of the gl's
%  pupil: defines pupil relative size (w.r.t. sSize), like a percentage
%  sSize: Size of cartesian coordinates. Space Size
%  disp_wrap: original (0) or wrapped mask (1)
%  plot_z: plot (1); no plot (0)
%  *-Not used directly here (but needed as an input of an inner function):
%  binMask: binarizes the mask w.r.t the max and min of the phase (boolean)
%  abs_ang: Magnitude (1); Phase (2) 
%  *-
%  monitorSize: size of the selected screen 
%  showM: show the mask. yes(1); no(0)
%
% Outputs:
%  Zernike phase mask
%
% Missing:
%  it still doesn't modulate on 2pi completely

%% Parameters
o1 = size(z_coeff); % Dimensions of z_coeff
o2 = abs(z_coeff) > ones(size(z_coeff)); % Ask if elemts are bigger than one
o2 = sum(o2); % Sum all elements; zero if they are all smaller than one
if z_coeff(1) == -1 && isscalar(z_coeff) == 1 % A); when z_coeff = -1
 z_vec = 2*pi*frac*a/L*randn(1,15)'; % Transposed random vector of integers
elseif o1(1) == 1 &&  sum(mod(z_coeff,1))== 0 % [#s] and integers % C)
 z_vec = zeros(1,15)'; % Transposed vector
 z_coeff =  z_coeff + 1; % Correction in order to get correctly the
                         % polynomials (j-th)
 z_vec(z_coeff) = 2*pi*frac*a/L; % Zernike weight vect: characterizes
                                 % the aberrations to be plotted
elseif o1(1) == 1 &&  o2 == 0 % Column vector and #s less than 1 [#s]' % B)
 z_vec = zeros(1,15)'; % Transposed vector
 for i = 1:size(z_coeff,2)
      z_vec(i) = 2*pi*frac*a*(z_coeff(i))/L;
 end
else
  warning('Not valid input');
  z_vec = zeros(1,2)'; % Initialization
end

%% Zernike Polynomials
n_mask = f_Zernike_Builder(z_vec,pupil,sSize,plot_z); % Defocus 
% (vector, pupil size, Matrix size (zernike phase size), graph:1 or not:0)
mask = exp(1i*n_mask); % Wrapped mask

%% Plot 
if disp_wrap == 1 && showM == 1
  figure; imagesc(x,y,angle(mask)), title('Wrapped phase mask');
  colormap(gray(gl)); cbh = colorbar; cbh.Label.String = 'Value of phase';
elseif showM == 1 % disp_wrap = 0
    abs_ang = 1; % "Magnitude" in order to not wrap the phase
    plotMask = showM; % plotMask = show; for 0 and 1.
    f_fig_maskSLM(x,y,r,n_mask,gl,glphi,mingl,maxgl,levShft,abs_ang, ...
                  binMask,monitorSize,plotMask);
    title('Unwrapped phase mask'); % The amplitude title is replaced
end

%% Phase wrapping test (optional): 
% The wrapped phase should do complete cycles of 2pi
% min(b,[],'omitnan') % Negative
% min(cos(b(:)),[],'omitnan') % positive since the phase was wrapped 
                              % for intensity   
end                              