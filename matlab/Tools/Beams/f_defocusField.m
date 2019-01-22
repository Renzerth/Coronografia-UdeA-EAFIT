function [defocusedImages,defocusedField,defocusMasks] = defocusField(complexField,lambda,observDistance,samplePlanes)
%% Computation Window definition
spaceSize = size(complexField);
xSLM = 2*spaceSize(1,1);
ySLM = 2*spaceSize(1,2); % computation window size 
[x2,y2]= meshgrid(-xSLM/2:xSLM/2-1,-ySLM/2:ySLM/2-1);
[~, rho2] = cart2pol(x2,y2);
apertureRatio = spaceSize(1,1);
Pupil = double(rho2 <= apertureRatio);
Pupil = Pupil.*exp(1i*2*pi*0);
%% Padding Window
coordX = ySLM/4+1:ySLM*3/4;
coordY = xSLM/4+1:xSLM*3/4;
u0_double = zeros(ySLM,xSLM);
u0_double(coordY,coordX) = complexField;
%% Defocus planes initialization
defocusMasks = zeros(ySLM,xSLM,samplePlanes);
PSF = zeros(ySLM,xSLM,samplePlanes);
defocusedField = zeros(ySLM/2,xSLM/2,samplePlanes);
defocusedImages = zeros(ySLM/2,xSLM/2,samplePlanes);
%% Defocus mask generation
delta = 3.75e-6;                        % pixel size for sensor, in m
delta_v = 8e-6;                         % pixel size for SLM, in m
k = 2*pi/lambda;
f = delta*delta_v*spaceSize(1,2)/lambda;       % In-focus distance m 
r_rec1 = sqrt(1-(lambda/delta*x2/xSLM).^2-(lambda/delta*y2/ySLM).^2); %% IMITATION OF VARIOUS DISTANCES IN SPACE DOMAIN
%% Defocused System
for index = 1:samplePlanes
  defocus = exp(1i*k*(observDistance*index*f)*r_rec1); % For large objects: observDistance*(index-1)+f
  PSF(:,:,index) = ifftshift(ifft2(fftshift(Pupil.*defocus)));
  defocusMasks(:,:,index) = fftshift(defocus);
end
%% Speckle defocused Images
speckleSpectrum = ifftshift(fft2(fftshift(u0_double)));
for index = 1:samplePlanes
  convKernelTF = ifftshift(fft2(fftshift(PSF(:,:,index))));
  u1 = ifftshift(ifft2(fftshift(convKernelTF.*speckleSpectrum)));
  defocusedImages(:,:,index) = abs(u1(coordY,coordX)).^2;
  defocusedField(:,:,index) = u1(coordY,coordX);
end
%% Plot generated intensities
fig_handler = figure(2);
for index = 1:samplePlanes
  figure(fig_handler);
  imagesc(defocusedImages(:,:,index)); axis square; colorbar;
  pause(0.20);
end
end