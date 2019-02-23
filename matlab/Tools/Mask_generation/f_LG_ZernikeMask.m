%% Laguerre Gauss Binary masks + Zernike Phases

function mask = f_LG_Zernike_Mask(x,y,r,phi,gl,glphi,mingl,maxgl, ...
                                  levShft,tc,s,ph0,p,W,binv,norm, ...
                                  abs_ang,z_coeff,a,frac,L,pupil,sSize, ...
                                  disp_wrap,plot_z,binMask,monitorSize, ...
                                  showM)
% Generates and plots a Laguerre Gauss + Zernike mask
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
%  p: number of radial nodes. If p=0, normal helicoid masks are obtained.
%     If they are used and tc=0(m=0); binary masks are obtained.
%     Even p; rings are ones. Odd p; rings are zeroes. Use mask = mask'
%  W: width of the modes: related with the radius of the phase and with the
%     disks on the magnitude
%  bininv: Binary inversion of the mask. Only applied when tc=0. Binary 
%          masks are only abtained when tc=0. yes(1); no(0)
%  norm: normalize magnitude and phase. yes(1); no(0)
%  abs_ang: Magnitude (1); Phase (2)
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
%  pupil: defines pupil relative size (w.r.t. sSize), like a percentage
%  sSize: Size of cartesian coordinates. Space Size
%  disp_wrap: original (0) or wrapped mask (1)
%  plot_z: plot (1); no plot (0)
%  binMask: binarizes the mask w.r.t the max and min of the phase (boolean)
%  monitorSize: size of the selected screen 
%  showM: show the mask. yes(1); no(0)
%
% Outputs:
% mask: VPL phase mask

  maskZ = f_Zernike_Mask(x,y,r,z_coeff,a,frac,L,gl,glphi,mingl,maxgl, ...
           levShft,pupil,sSize,disp_wrap,plot_z,binMask,monitorSize,showM);
  maskLG = f_LG_Mask(x,y,r,phi,gl,glphi,mingl,maxgl,levShft,tc,s,ph0,p, ...
                     W,binv,norm,abs_ang,binMask,monitorSize,showM);
  mask = maskZ.*maskLG; % Combined mask. It is recommended to use binary
                        % masks on LG
  
  %% Plot the combined LG + Z mask                      
  wrappedMask = f_mask_circ_angle_gl(r,mask,binMask,glphi,mingl,maxgl, ...
                                     levShft);
  tit = strcat('LG phase mask with topological charge ', ...
         num2str(tc),' and radial node ',num2str(p),' + Zernike mask');
  f_fig_maskPCscreen(x, y, wrappedMask, tit, gl, showM);  
end