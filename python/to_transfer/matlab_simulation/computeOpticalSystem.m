function [OTF,PSF,MTF,phaseTF] = computeOpticalSystem(systemPupil,dataCoordinates,analysisScaling,synthesisScaling)
%% Imaging System Pupil properties
systemResponse = ifftshift(ifft2(fftshift(systemPupil)))*synthesisScaling;
windowedOTF = ifftshift(fft2(fftshift(conj(systemResponse).*systemResponse)))*analysisScaling; % Incoherent transfer function -- Selfcorrelation = selfconvolution of pupil
OTF = windowedOTF(dataCoordinates(1,:),dataCoordinates(2,:));
%% Spectral response properties
impulseResponse = ifftshift(ifft2(fftshift(windowedOTF)))*synthesisScaling; % Incoherent impulse response
PSF = abs(impulseResponse(dataCoordinates(1,:),dataCoordinates(2,:))).^2;
MTF = abs(OTF)./(max(abs(OTF(:))));
phaseTF = angle(OTF)./(max(angle(OTF(:))));
end