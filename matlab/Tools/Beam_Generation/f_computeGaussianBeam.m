function [gaussianBeam,x,y,k] = computeGaussianBeam(sizeRatio,spaceSize,beamCenter,lambda)
%% Space definition
xLength = spaceSize(1,1);
yLength = spaceSize(1,2);
xCoordinate = -xLength/2:xLength/2-1;
yCoordinate = -yLength/2:yLength/2-1;
[x,y]= meshgrid(xCoordinate,yCoordinate);
%% Gaussian Parameters % e.g. Orazio Svelto, Principles of Lasers, 4th ed. page 152, Eqn's 4.715a-4.717c / OSLO Ref
w0 = sizeRatio*xLength/2;      %Beam field width at waist
z = 1e-20;                      %Distance from beam axis
xC = beamCenter(1,1);           %x Gaussian center
yC = beamCenter(1,2);           %y Gaussian center

k  = 2*pi/lambda;              %Wavenumber
zR = pi*w0^2/lambda;           %Rayleigh range
b0 = 2*zR;                     %Confocal Radius

w = w0*sqrt(1+(z/zR)^2);      %Beam (field) width
R = z*(1+(zR/z)^2);           %Radius of curvature
phi = atan(z/zR);             %Guoy phase

gaussianBeam = (w0/w) * exp(-((x - xC).^2 + (y - yC).^2)/w^2) .* exp(-1i*k*((x - xC).^2 + (y - yC).^2)/2/R) * exp(1i*phi);
end