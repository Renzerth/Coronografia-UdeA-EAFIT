% Angular Spectrum approximation allows both sampling rate and space size
% to be equal. It remains valid as long the pattern's size does not vary
% significantly from the near-field.
%
% Space quadratic phase factors-based techniques require a multisampling
% rate due to its higher frequencies sampling which can be satisfied through
% Nyquist sampling in the space scaling factors thus space size is modified.
% If it is held constant, bad sampling will lead to aliased Fresnel intengral.
% A padding of the optical field will increase its sampling though.
% However, at long distance field's  phase high variation will limit solution
% to intensity-only analysis due to large undersampling of phase.
%
% Further reading: Optics Communications 164 (1999). 233â€?45.

%% Spectral sampling settings
samplingFactor = 1;
enablePadding = samplingFactor -1;
%% Plane properties
planeSize = [0.5,0.5]; % m
propDistance = 100; % m
%% Light properties
lambda = 0.6328e-6; % Light central wavelength [m]
coherentLight = true;
%% System Properties
TC = 2;
pupilRadius = 1;
kernelType = 1;
focalLength = 1000; % [m]
% aberrationIndex = 4;
%% Gaussian beam parameters
gaussianSize = 0.1;
spaceSize = [500, 500];
halfSize = ceil((spaceSize+1)/2);
beamCenter = [0,0];
%% Uniform light parameters
circSize = 0.5;
%% Beam coordinates
[gaussianBeam, X, Y, k] = BEAMS.computeGaussianBeam(gaussianSize,spaceSize,beamCenter,lambda);
rho = sqrt( X.^2 + Y.^2);
phi = atan2(Y, X);
circBeam = double(rho <= halfSize(1)*circSize);
%% Spectral properties
[freqVectX, freqVectY,spXperiod,spYperiod,analysisScaling,normNMFactor,synthesisScaling] = computeFreqVector(planeSize, spaceSize, samplingFactor);
%% Input FIeld
inputPlane = gaussianBeam;
%% Lens properties
lensRadius = sqrt((spXperiod*X).^2 + (spYperiod*Y).^2);
f = spXperiod*spYperiod*max(spaceSize)/lambda;       % In-focus distance m
lensPhaseA = exp(1i*k*focalLength*(sqrt(1-(lambda/spXperiod*X/spaceSize(1)).^2-(lambda/spYperiod*Y/spaceSize(2)).^2)));
lensPhaseB = exp(-1i*k/(2*focalLength)*(rho/spaceSize(1)).^2);
lensPhaseC = exp(-1i*k/(2*focalLength)*(lensRadius).^2);
%% Prepare propagator
[beamSpectrum, dataCoordX, dataCoordY, freqMeshX, freqMeshY,freqCirc] = prepareSpectralProp(spaceSize,samplingFactor, inputPlane, freqVectX, freqVectY, analysisScaling, lambda);
[propagationKernel] = computeConvolutionKernel(freqMeshX,freqMeshY,freqCirc,k, propDistance, planeSize,kernelType);
%% Optical System ZernikePhase
% zernikeCoeffs = zeros(1,15);
% zernikeCoeffs(aberrationIndex) = 0;
% [systemPhase,~] = ZERNIKE.Zernike_Builder(zernikeCoeffs',pupilRadius,min(spaceSize),false);
% systemPhase(isnan(systemPhase)) = 0;
Pupil = double(rho <= halfSize(1)*pupilRadius);
%% System Definition
systemPupil = zeros(samplingFactor*spaceSize);
systemPupil(dataCoordX,dataCoordY) = Pupil.*exp(1i*TC*phi).*lensPhaseC;
[OTF,~,~,~] = computeOpticalSystem(systemPupil, [dataCoordX; dataCoordY], analysisScaling,synthesisScaling);
%% Coherence Response Definition
if coherentLight
    systemTransferFunct = fftshift(systemPupil);
else
    systemTransferFunct = fftshift(padarray(OTF,halfSize*enablePadding,0,'both'));
end
%% Propagation
propagationKernel = fftshift(propagationKernel);
spectralProduct = beamSpectrum.*propagationKernel.*systemTransferFunct;
filteredSignal = ifft2s(spectralProduct,normNMFactor,synthesisScaling);
inRangeData = filteredSignal(dataCoordX,dataCoordX);
%% Field Properties
outputIntensity = abs(inRangeData).^2;
inputIntensity = abs(inputPlane).^2;
%% PLots
close all;

figure('Color', 'white');
subplot(2,3,1); imagesc(inputIntensity); title('Input Plane'); axis square;

subplot(2,3,2); imagesc(outputIntensity); title('Output Plane'); axis square;

subplot(2,3,3); imagesc(angle(ifftshift(systemTransferFunct))); title('Optical System TF Phase'); axis square;

subplot(2,3,4); imagesc(log10(inputIntensity)); title('Input Log10 view'); axis square;

subplot(2,3,5); imagesc(log10(outputIntensity)); title('Output Log10 view'); axis square;

subplot(2,3,6); imagesc(angle(ifftshift(propagationKernel))); title('Propagation Transfer Function'); axis square;