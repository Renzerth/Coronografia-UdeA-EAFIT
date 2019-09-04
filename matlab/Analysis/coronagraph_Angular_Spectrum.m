%% Plane properties
planeSize = 2*[2.54e-2, 2.54e-2]; % m
spaceSamples = [1024, 1024];
halfSize = floor((spaceSamples+1)/2);
shiftDistance = 0;

%% Light properties
lambda = 0.6328e-6; % Light central wavelength [m]
k = 2*pi/lambda;

insidentEnergy = 1; % [W/m^2]
mediumElecPermitivity = 8.85e-12; % permittivity of free space [F/m]
mediumRefracIndex = 1; % refractive index
propagationSpeed = 3e8; % speed of light [m/s]
irradianceScaling = (propagationSpeed*mediumRefracIndex*mediumElecPermitivity)/2; % [W/m^2]

%% System Properties
illuminationDiameter = 2e-3; %[m]
LyotApertureDiameter = illuminationDiameter; % [m]
aberrationPupilRadius = max(planeSize); %[m]

focalLengthA = 2; % [m]
focalLengthB = 2;
aberrationIndex = 8;
aberrationValue = 0*pi*(1.0021)*0.01;

% FresnelNumber = (2*pupilRadius)^2/(propDistance*focalLengthA);

%% Spatial coordinates
xCenter= floor((spaceSamples(1) + 1)/2);
yCenter = floor((spaceSamples(2) + 1)/2);
[xCoords, yCoords, X, Y, theta, rho] = makeSpaceCoordinates(planeSize(1),planeSize(2),spaceSamples(1),spaceSamples(2),xCenter,yCenter);

%% Spectral coordinates
samplingFactor = 1;
[freqVectX, freqVectY,spXperiod,spYperiod,analysisScaling,normNMFactor,synthesisScaling] = computeFreqVector(planeSize, spaceSamples, samplingFactor);

%% Vortex Mask properties
TC = 0;
grayLevels = 256;
phaseValues = linspace(-pi,pi,grayLevels);
% vortexMask = TOOLS.spiralGen2(spaceSamples,TC);
vortexMask = exp(1i*TC*theta);
[discretMap,~] = TOOLS.discretizeMask(phaseValues,angle(vortexMask));
% vortexMask = exp(1i*discretMap);

%% Lens properties
lensRadii = rho;
lensDiameter = 2.54e-2;
lensAperture = double(rho <= lensDiameter/2);
lensA = lensAperture.*exp(-1i*k/(2*focalLengthA)*(lensRadii).^2);
lensB = lensAperture.*exp(-1i*k/(2*focalLengthB)*(lensRadii).^2);
diverLens = lensAperture.*exp(1i*k/(2*focalLengthA)*(lensRadii).^2);

%% Uniform light definition after first focal lens
inputPlane = insidentEnergy*double(rho <= illuminationDiameter/2).*lensA;

%% Optical Aberrations Zernike Phase
zernikeCoeffs = zeros(1,15);
zernikeCoeffs(aberrationIndex) = aberrationValue;
[systemPhase,~] = Zernike_Builder(zernikeCoeffs',aberrationPupilRadius/planeSize(1),min(spaceSamples),false);
systemPhase(isnan(systemPhase)) = 0;

%% Optical System Phase Properties Definition
Pupil = double(rho <= aberrationPupilRadius);
pupilTransferFunct = Pupil.*exp(1i*systemPhase);

%% Prepare propagators
propDistance = focalLengthA - shiftDistance;
[dataCoordX, dataCoordY, freqMeshX, freqMeshY, freqRadii, freqCirc] = prepareSpectralProp_mod(spaceSamples, freqVectX, freqVectY, lambda);
distanceShiftA = exp(-1i*k*focalLengthA);
distanceShiftB = exp(-1i*k*focalLengthB);

[focalPlanePropagationKernelA] = freqCirc.*exp(1i*k*focalLengthA*sqrt(1 - (lambda*freqMeshX).^2 - (lambda*freqMeshY).^2));
[focalPlanePropagationKernelB] = freqCirc.*exp(1i*k*focalLengthB*sqrt(1 - (lambda*freqMeshX).^2 - (lambda*freqMeshY).^2));
[shiftDistanceKernel] = freqCirc.*exp(1i*k*propDistance*sqrt(1 - (lambda*freqMeshX).^2 - (lambda*freqMeshY).^2));

%% Generate shifted focal distribution - SLM Input
systemTransFunct = pupilTransferFunct.*shiftDistanceKernel;
[focalPlaneDistribution] = distanceShiftA*convoluteSignal(systemTransFunct, inputPlane, normNMFactor, analysisScaling, synthesisScaling, dataCoordX, dataCoordY, halfSize);

%% Vortex Focal Plane Filtering
SLMfilteredField = focalPlaneDistribution.*vortexMask;

%% Generate Lyot distribution - Pupil Plane
[ocularPlaneDistribution] = distanceShiftA*deconvoluteSignal(focalPlanePropagationKernelA, SLMfilteredField, normNMFactor, analysisScaling, synthesisScaling, dataCoordX, dataCoordY, halfSize);
LyotPlaneDistribution = ocularPlaneDistribution.*diverLens;

%% Lyot Pupil Truncation - Pupil Plane
LyotAperture = double(rho < LyotApertureDiameter/2);
truncatedLyotPlane = distanceShiftA*LyotPlaneDistribution.*LyotAperture;

%% Add Analytic Lens Converging Phase
PSFlensPlaneDistribution = distanceShiftB*truncatedLyotPlane.*lensB;

%% Generate PSF response focal point
[PSFplaneDistribution] = distanceShiftB*deconvoluteSignal(focalPlanePropagationKernelB, PSFlensPlaneDistribution, normNMFactor, analysisScaling, synthesisScaling, dataCoordX, dataCoordY, halfSize);

%% Field Properties
compIntensity = @(complexField, energyScaling) energyScaling*abs(complexField).^2;
compAngle = @(complexField) angle(complexField);

outputIntensity = compIntensity(PSFplaneDistribution,irradianceScaling);
inputIntensity = compIntensity(inputPlane,insidentEnergy);

%% PLots
close all;
apertureDisplayRatio = 1.5;
PSFdisplayRatio = 0.25;
dispRngA = {abs(xCoords) <= illuminationDiameter*apertureDisplayRatio; abs(yCoords) <= illuminationDiameter*apertureDisplayRatio};
dispRngB =  {abs(xCoords) <= planeSize(1)*PSFdisplayRatio; abs(yCoords) <= planeSize(2)*PSFdisplayRatio};

figureHandleA = figure('Color', 'white');
subplot(2,3,1); imagesc(xCoords,yCoords,inputIntensity); title('Input Plane'); xlabel('[m]'); ylabel('[m]'); CH = colorbar; CH.Label.String = 'W/m^{2}'; axis square;
subplot(2,3,2); imagesc(xCoords,yCoords,outputIntensity); title('Output Plane'); xlabel('[m]'); ylabel('[m]'); CH = colorbar; CH.Label.String = 'W/m^{2}'; axis square;
subplot(2,3,3); imagesc(xCoords,yCoords,angle(inputPlane)); title(sprintf('L1 Lens Phase: f=%1.1fm',focalLengthA)); xlabel('[m]'); ylabel('[m]'); CH = colorbar; CH.Label.String = 'rad'; axis square;
subplot(2,3,4); imagesc(xCoords,yCoords,10*log10(inputIntensity)); title('Input Log10 view'); xlabel('[m]'); ylabel('[m]'); CH = colorbar; CH.Label.String = 'dB'; axis square;
subplot(2,3,5); imagesc(xCoords,yCoords,10*log10(outputIntensity)); title('Output Log10 view'); xlabel('[m]'); ylabel('[m]');CH = colorbar; CH.Label.String = 'dB'; axis square;
subplot(2,3,6); imagesc(freqVectX,freqVectY,angle(focalPlanePropagationKernelA)); title(sprintf('Vaccum Transfer Function - D=%1.1fm', propDistance)); xlabel('fx [1/m]'); ylabel('fy [1/m]'); CH = colorbar; CH.Label.String = 'rad'; axis square;
set(gcf,'units','normalized','position',[0,0,1,1]);
% TOOLS.tightfigadv(figureHandle);

figureHandleB = figure('Color', 'white');
subplot(2,3,1); imagesc(xCoords(dispRngA{1}),yCoords(dispRngA{2}),inputIntensity(dispRngA{1},dispRngA{2})); title('Input Plane'); xlabel('[m]'); ylabel('[m]'); CH = colorbar; CH.Label.String = 'W/m^{2}'; axis square;
subplot(2,3,2); imagesc(xCoords(dispRngA{1}),yCoords(dispRngA{2}),compIntensity(focalPlaneDistribution(dispRngA{1},dispRngA{2}),irradianceScaling)); title(sprintf('Focal Plane after L1 - EFL:%1.1fm.',focalLengthA)); xlabel('[m]'); ylabel('[m]'); CH = colorbar; CH.Label.String = 'W/m^{2}'; axis square;
subplot(2,3,3); imagesc(xCoords(dispRngA{1}),yCoords(dispRngA{2}),compIntensity(ocularPlaneDistribution(dispRngA{1},dispRngA{2}),irradianceScaling)); title(sprintf('Plane at L2 at entrance - z = %1.1fm',2*focalLengthA)); xlabel('[m]'); ylabel('[m]'); CH = colorbar; CH.Label.String = 'W/m^{2}'; axis square;
subplot(2,3,4); imagesc(xCoords(dispRngA{1}),yCoords(dispRngA{2}),compIntensity(truncatedLyotPlane(dispRngA{1},dispRngA{2}),irradianceScaling)); title(sprintf('Truncated Lyot Distribution - %1.1fmm.', LyotApertureDiameter*10^abs(ceil(log10(LyotApertureDiameter)-1)))); xlabel('[m]'); ylabel('[m]'); CH = colorbar; CH.Label.String = 'W/m^{2}'; axis square;
subplot(2,3,5); imagesc(xCoords(dispRngA{1}),yCoords(dispRngA{2}),compIntensity(PSFlensPlaneDistribution(dispRngA{1},dispRngA{2}),irradianceScaling)); title(sprintf('Plane at L3 entrance - EFL:%1.1fm.',focalLengthB)); xlabel('[m]'); ylabel('[m]'); CH = colorbar; CH.Label.String = 'W/m^{2}'; axis square;
subplot(2,3,6); imagesc(xCoords(dispRngB{1}),yCoords(dispRngB{2}),outputIntensity(dispRngB{1},dispRngB{2})); title(sprintf('PSF Distribution - z=%1.1fm.',3*focalLengthA + 2*focalLengthB)); xlabel('[m]'); ylabel('[m]'); CH = colorbar; CH.Label.String = 'W/m^{2}'; axis square;
set(gcf,'units','normalized','position',[0,0,1,1]);

figureHandleC = figure('Color', 'white');
subplot(2,3,1); imagesc(xCoords,yCoords,compAngle(inputPlane)); title('Input Plane'); xlabel('[m]'); ylabel('[m]'); CH = colorbar; CH.Label.String = 'rad'; axis square;
subplot(2,3,2); imagesc(xCoords,yCoords,compAngle(focalPlaneDistribution)); title(sprintf('Focal Plane after L1 - EFL:%1.1fm.',focalLengthA)); xlabel('[m]'); ylabel('[m]'); CH = colorbar; CH.Label.String = 'rad'; axis square;
subplot(2,3,3); imagesc(xCoords,yCoords,compAngle(ocularPlaneDistribution)); title(sprintf('Plane at L2 at entrance - z = %1.1fm.',2*focalLengthA)); xlabel('[m]'); ylabel('[m]'); CH = colorbar; CH.Label.String = 'rad'; axis square;
subplot(2,3,4); imagesc(xCoords,yCoords,compAngle(truncatedLyotPlane)); title(sprintf('Truncated Lyot Distribution - %1.1fmm.', LyotApertureDiameter*10^abs(ceil(log10(LyotApertureDiameter)-1)))); xlabel('[m]'); ylabel('[m]'); CH = colorbar; CH.Label.String = 'rad'; axis square;
subplot(2,3,5); imagesc(xCoords,yCoords,compAngle(PSFlensPlaneDistribution)); title(sprintf('Plane at L3 entrance - EFL:%1.1fm.',focalLengthB)); xlabel('[m]'); ylabel('[m]'); CH = colorbar; CH.Label.String = 'rad'; axis square;
subplot(2,3,6); imagesc(xCoords,yCoords,compAngle(PSFplaneDistribution)); title(sprintf('PSF Distribution - z=%1.1fm.',3*focalLengthA + 2*focalLengthB)); xlabel('[m]'); ylabel('[m]'); CH = colorbar; CH.Label.String = 'rad'; axis square;
set(gcf,'units','normalized','position',[0,0,1,1]);