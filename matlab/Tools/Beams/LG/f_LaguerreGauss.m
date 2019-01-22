function [ LG ] = f_LaguerreGauss(r, phi, m, s, ph0, p, W)
% On axis plot of the Laguerre-Gauss or Optical Vortex Beam
%
% Inputs
%  rho: radius in cylindrical coordinates
%  phi: angle in cylindrical coordiantes
%  m: order of the derivative that defines it; positive for the radial part
%     but the sign can vary for the azimuthal: direction of the spiral 
%     mask. Correspons to the topological charge of the spiral phase mask
%     or the phase discontinuity. It is known as the azimuthal index giving
%     an OAM of l*planck_reduced per photon. Vorticity increases with it.
%  s: sign of spiral phase mask (+1 or -1)
%  ph0: initial phase of the azimuthal part (and of the spiral phase mask)
%  p: degree of polynomial. Is known as the number of radial nodes in the
%     intensity distribution. There are p+1 rings or airy disks present. 
%     Denotes the rings transitions
%  w = w(0): radius of the beam; could depend on z: decreases along the
%            optical axis for example
%
% Notes:
%  z = 0; an incidence plane is assumed during calculations
%
% LG beams taken from: 
%  1_book_Orbital angular momentum origins_behavior_applications_2011
%  Page: 169; it has z dependence (not the case here)

normFactor = sqrt( 2*factorial(p) / (pi*W.^2*factorial(p+abs(m))) );
rNorm = r*sqrt(2)*W.^(-1);
U = exp(-rNorm.^2/2); % Gaussian function
radial = (rNorm).^abs(m) .* U .* ...
 polyval(f_LaguerreAssociated(p, abs(m)), rNorm.^2); % evaluates the pol.
azimuthal = exp( 1i*m*s*(phi+ph0) );
LG = normFactor .* radial .* azimuthal; % Output as a complex field
end