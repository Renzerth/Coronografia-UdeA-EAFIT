function angleAiry = f_scalePix2DiffAng(x,AiryFactor)
% Cconverts from pixels to the diffraction angle in Lamda/D units
%
% This function only converts pixels to L/D. Applies only on one Airy disk
%
% R = 1.22*L/(2*NA) = 1.22*L*f/(n*D) = x*PP*M -> L/D = (n*PP*M)/(1.22*f)*x, 
% with m*sin(theta) ~ m*theta =  m*L/D [m=1 is the first diffraction order]
% 
% Inputs:
%  AiryFactor = n*PP*M/f; The calculation were explained above [Lambda/D]
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
BesselFirstZero = 1.22; % First zero of the cylindrical Bessel function of 
                        % first kind and zeroth order
angleAiry = (AiryFactor/BesselFirstZero)*x; % [Lambda/D] Units
end

