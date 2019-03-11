function angleAiry = f_scalePix2DiffAng(x,n,PP,M,f)
% Inputs:
%  x: vector of a pixel distance that comes from a camera image
%  n: refractive index of the medium
%  PP: pixel pitch of the camera [um]
%  M: magnification of the camera's objective
%  f: focal distance of the camera's objective [mm]
%  
% Output:
%  angleAiry: L/D (wavelength over system's diameter), angular units in
%  rad. A vector scaled from pixel size to diffraction angle size
%
% R = 1.22*L/(2*NA) = 1.22*L*f/(n*D) = x*PP*M -> L/D = (n*PP*M)/(1.22*f)*x, 
% with m*sin(theta) ~ m*theta =  m*L/D [m=1 is the first diffraction order]
f = f*1e3; % mm to um
angleAiry = (n*PP*M)/(1.22*f)*x;
end

