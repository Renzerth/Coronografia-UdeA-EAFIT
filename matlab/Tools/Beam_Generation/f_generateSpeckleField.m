function [speckleImage,U0] = f_generateSpeckleField(sizeRatio,spaceSize,beamCenter,lambda,roughnessSTD,roughnessMean,pupilRatio)
%% Gaussian Beam 
[gaussianBeam,x,y,k] = f_computeGaussianBeam(sizeRatio,spaceSize,beamCenter,lambda);
gaussianIntensity = abs(gaussianBeam).^2;
[~, rho] = cart2pol(x,y);
%% Rough Surface definition
roughnessDistb = roughnessSTD*randn(spaceSize(1,1),spaceSize(1,2)) + roughnessMean; %rng('shuffle')
surfacePhaseFactor = exp(1i*2*k*roughnessDistb);
u0 = gaussianBeam.*surfacePhaseFactor;
%% Imaging System Pupil properties
apertureRatio = pupilRatio*spaceSize(1,1)/2;
Pupil = double(rho<=apertureRatio);
Pupil = Pupil.*exp(1i*2*pi*0); % Coherent transfer function
cohImpulseResponse = ifftshift(ifft2(fftshift(Pupil)));
%% Speckle Image by convolution
fieldSpectrum = ifftshift(fft2(fftshift(u0)));
convKernelTF = ifftshift(fft2(fftshift(cohImpulseResponse)));
U0 = ifftshift(ifft2(fftshift(convKernelTF.*fieldSpectrum)));
speckleImage = abs(U0).^2;
%% Plot
figure('name','Speckle generation by imaging rough surface','units','normalized','outerposition',[0 0 1 1],'color','white');

subplot(1,3,1);
imagesc(gaussianIntensity); axis square;
titleTextA = sprintf('Gaussian intensity distribution at: %i%% of width.',100*sizeRatio);
title(titleTextA,'FontSize',15,'FontWeight','bold');
xlabel('X pixel number','FontSize',15,'FontWeight','normal')
ylabel('Y pixel number','FontSize',15,'FontWeight','normal')
axis square;

subplot(1,3,2);
imagesc(speckleImage); axis square;
titleTextB = sprintf('Speckle intensity at: %i%% of pupil size.',100*pupilRatio);
title(titleTextB,'FontSize',15,'FontWeight','bold');
xlabel('X pixel number','FontSize',15,'FontWeight','normal')
ylabel('Y pixel number','FontSize',15,'FontWeight','normal')
axis square;

subplot(1,3,3);
imagesc(angle(U0)); axis square;
title('Speckle field phase','FontSize',15,'FontWeight','bold');
xlabel('X pixel number','FontSize',15,'FontWeight','normal')
ylabel('Y pixel number','FontSize',15,'FontWeight','normal')
axis square;
end
