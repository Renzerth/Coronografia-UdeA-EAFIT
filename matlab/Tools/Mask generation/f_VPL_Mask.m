%% Vortex-Producing Lens (VPL) phase mask
% Taken from
% 1_edgar_2013_High-quality optical vortex-beam generation_E-Rueda_OL.pdf
% Equation 3, page 2

function mask = f_VPL_Mask(x,y,r,phi,gl,tc,s,ph0,L,f_FR,binMask,showM)
% Generates and plots a VPL mask:  helicoidal mask + fresnel lens
%
% Inputs: 
%  x,y: cartesian coordinates vector
%  r,phi: polar coordinates (r in cm)
%  gl: number of grey levels (normally 256)
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
% mask: Vortex-Producing Lens (VPL) phase mask

    %% Spiral phase mask Generation
    m = s*tc;  % tc with a sign
    maskSPP = m*(phi + ph0); % General mask. Angle phi is wrapped on [-pi,pi]
    maskSPP = exp(1i*maskSPP); % Wrapped mask

    %% VPL Mask
    r = r*1e4; % r is a in cm and converted to um
    VPLfactor = -(pi*r.^2)/(L*f_FR);
    maskVPL = exp(1i*VPLfactor);
    mask = maskSPP.*maskVPL; % mask*maskVPL; gives cool results as the Fresnel
                          % lens but in a lineal fashion: a chirp ramp for
                          % tc = 0

%% Plot (with axes)
if showM == 1
  wrappedMask = f_circularPupil_maskAngle(r,mask,binMask);
  tit = ['VPL with topological charge ' ...
         num2str(tc) ' and ' num2str(gl) ' gray levels'];  
  f_fig_maskPCscreen(x, y, wrappedMask, tit, gl, showM);
end

end