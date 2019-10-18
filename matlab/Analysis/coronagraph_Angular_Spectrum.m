%% path dependencies
addpath(genpath(fileparts(pwd)));

%% Plane properties
planeSize = 2*[2.54e-2, 2.54e-2]; % m
spaceSamples = 2*[1024, 1024];
halfSize = floor((spaceSamples+1)/2);
shiftDistance = 0;

%% Light properties
lambda = 0.6328e-6; % Light central wavelength [m]
k = 2*pi/lambda;

insidentEnergy = 1; % [W/m^2]
mediumElecPermitivity = 8.85e-12; % permittivity of free space [F/m]
mediumRefracIndex = 1; % refractive index
propagationSpeed = 3e8; % speed of light [m/s]
irradianceScaling = insidentEnergy*(propagationSpeed*mediumRefracIndex*mediumElecPermitivity)/2; % [W/m^2]

%% System Properties
illuminationDiameter = 6e-3; %[m]
obstructionDiameter = illuminationDiameter/2;
LyotApertureDiameter = illuminationDiameter; % [m]
aberrationPupilRadius = max(planeSize); %[m]

focalLengthA = 2; % [m]
focalLengthB = 2;
aberrationIndex = 8;
aberrationValue = 0*pi*(1.0021)*0.01;
FresnelNumber = (illuminationDiameter)^2/(focalLengthA*lambda);

%% Spatial coordinates
xCenter= floor((spaceSamples(1) + 1)/2);
yCenter = floor((spaceSamples(2) + 1)/2);
[xCoords, yCoords, X, Y, theta, rho] = makeSpaceCoordinates(planeSize(1),planeSize(2),spaceSamples(1),spaceSamples(2),xCenter,yCenter);

%% Spectral coordinates
samplingFactor = 1;
[freqVectX,freqVectY,spXperiod,spYperiod,analysisScaling,normNMFactor,synthesisScaling] = computeFreqVector(planeSize, spaceSamples, samplingFactor);

%% Vortex Mask properties
TC = 5;
grayLevels = 24;
vortexMask = spiralGen2(spaceSamples,TC);
% vortexMask = exp(1i*(TC*(theta)));
[vortexMask] = discretizeMap(angle(vortexMask),grayLevels);

%% Lens properties
lensRadii = rho;
lensDiameter = 2.54e-2;
lensAperture = double(rho <= lensDiameter/2);
lensA = lensAperture.*exp(-1i*k/(2*focalLengthA)*(lensRadii).^2);
lensB = lensAperture.*exp(-1i*k/(2*focalLengthB)*(lensRadii).^2);
diverLens = lensAperture.*exp(1i*k/(2*focalLengthA)*(lensRadii).^2);

%% Uniform light definition after first focal lens
% aperture = double(rho <= illuminationDiameter/2);
aperture = (rho <= illuminationDiameter/2).*(rho >= obstructionDiameter/2); % Obstructed
inputPlane =  insidentEnergy*aperture.*lensA;
% inputPlane = insidentEnergy*aperture.*evaluateGaussianField(X,Y,illuminationDiameter/3,[0,0],lambda).*lensA;

%% Optical Aberrations Zernike Phase
zernikeCoeffs = zeros(1,10);
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
[PSFplaneDistribution] = distanceShiftB*convoluteSignal(focalPlanePropagationKernelB, PSFlensPlaneDistribution, normNMFactor, analysisScaling, synthesisScaling, dataCoordX, dataCoordY, halfSize);

%% Field Properties
compIntensity = @(complexField, energyScaling) energyScaling*abs(complexField).^2;
compAngle = @(complexField) angle(complexField);

outputIntensity = compIntensity(PSFplaneDistribution,irradianceScaling);
inputIntensity = compIntensity(inputPlane,irradianceScaling);

%% PLots Settings
close all;
fontSize = 16;
apertureDisplayRatio = 1.5;
PSFdisplayRatio = 0.1;

plotGaps = [0.1, 0.01];
heighMargins = [0.1, 0.05];
widthMargins = [.02, 0.0];

dispRngA = {abs(xCoords) <= illuminationDiameter*apertureDisplayRatio; abs(yCoords) <= illuminationDiameter*apertureDisplayRatio};
dispRngB =  {abs(xCoords) <= planeSize(1)*PSFdisplayRatio; abs(yCoords) <= planeSize(2)*PSFdisplayRatio};

%% Plots
figureHandleA = figure('Color', 'white');
[ha, ~] = tight_subplot(2,3,plotGaps,heighMargins,widthMargins);
axes(ha(1)); imagesc(xCoords,yCoords,inputIntensity); title('Input Plane'); xlabel('[m]','FontSize',fontSize,'FontWeight','bold'); ylabel('[m]','FontSize',fontSize,'FontWeight','bold'); CH = colorbar; CH.Label.String = 'W/m^{2}'; axis square; set(gca,'FontSize',fontSize,'FontWeight','normal');
axes(ha(2)); imagesc(xCoords,yCoords,outputIntensity); title('Output Plane'); xlabel('[m]','FontSize',fontSize,'FontWeight','bold'); ylabel('[m]','FontSize',fontSize,'FontWeight','bold'); CH = colorbar; CH.Label.String = 'W/m^{2}'; axis square; set(gca,'FontSize',fontSize,'FontWeight','normal');
axes(ha(3)); imagesc(xCoords,yCoords,angle(inputPlane)); title(sprintf('L1 Lens Phase: f=%1.1fm',focalLengthA)); xlabel('[m]','FontSize',fontSize,'FontWeight','bold'); ylabel('[m]','FontSize',fontSize,'FontWeight','bold'); CH = colorbar; CH.Label.String = 'rad'; axis square; set(gca,'FontSize',fontSize,'FontWeight','normal');
axes(ha(4)); imagesc(xCoords,yCoords,10*log10(inputIntensity)); title('Input Log10 view'); xlabel('[m]','FontSize',fontSize,'FontWeight','bold'); ylabel('[m]','FontSize',fontSize,'FontWeight','bold'); CH = colorbar; CH.Label.String = 'dB'; axis square; set(gca,'FontSize',fontSize,'FontWeight','normal');
axes(ha(5)); imagesc(xCoords,yCoords,10*log10(outputIntensity)); title('Output Log10 view'); xlabel('[m]','FontSize',fontSize,'FontWeight','bold'); ylabel('[m]','FontSize',fontSize,'FontWeight','bold');CH = colorbar; CH.Label.String = 'dB'; axis square; set(gca,'FontSize',fontSize,'FontWeight','normal');
axes(ha(6)); imagesc(freqVectX,freqVectY,angle(focalPlanePropagationKernelA)); title(sprintf('Vaccum Transfer Function: D=%1.1fm', propDistance)); xlabel('fx [1/m]','FontSize',fontSize,'FontWeight','bold'); ylabel('fy [1/m]','FontSize',fontSize,'FontWeight','bold'); CH = colorbar; CH.Label.String = 'rad'; axis square; set(gca,'FontSize',fontSize,'FontWeight','normal');
set(figureHandleA,'units','normalized','position',[0,0,1,1]);

figureHandleB = figure('Color', 'white');
[ha, ~] = tight_subplot(2,3,plotGaps,heighMargins,widthMargins);
axes(ha(1)); imagesc(xCoords(dispRngA{1}),yCoords(dispRngA{2}),inputIntensity(dispRngA{1},dispRngA{2})); title('Input Plane','FontSize',fontSize,'FontWeight','bold'); xlabel('[m]','FontSize',fontSize,'FontWeight','bold'); ylabel('[m]','FontSize',fontSize,'FontWeight','bold'); CH = colorbar; CH.Label.String = 'W/m^{2}'; axis square; set(gca,'FontSize',fontSize,'FontWeight','normal');
axes(ha(2)); imagesc(xCoords(dispRngA{1}),yCoords(dispRngA{2}),compIntensity(focalPlaneDistribution(dispRngA{1},dispRngA{2}),irradianceScaling)); title(sprintf('Focal Plane: EFL=%1.1fm',focalLengthA)); xlabel('[m]','FontSize',fontSize,'FontWeight','bold'); ylabel('[m]','FontSize',fontSize,'FontWeight','bold'); CH = colorbar; CH.Label.String = 'W/m^{2}'; axis square; set(gca,'FontSize',fontSize,'FontWeight','normal');
axes(ha(3)); imagesc(xCoords(dispRngA{1}),yCoords(dispRngA{2}),compIntensity(ocularPlaneDistribution(dispRngA{1},dispRngA{2}),irradianceScaling)); title(sprintf('Plane at L2: z=%1.1fm',2*focalLengthA)); xlabel('[m]','FontSize',fontSize,'FontWeight','bold'); ylabel('[m]','FontSize',fontSize,'FontWeight','bold'); CH = colorbar; CH.Label.String = 'W/m^{2}'; axis square; set(gca,'FontSize',fontSize,'FontWeight','normal');
axes(ha(4)); imagesc(xCoords(dispRngA{1}),yCoords(dispRngA{2}),compIntensity(truncatedLyotPlane(dispRngA{1},dispRngA{2}),irradianceScaling)); title(sprintf('Truncated Lyot: %1.1fmm', LyotApertureDiameter*10^abs(ceil(log10(LyotApertureDiameter)-1)))); xlabel('[m]','FontSize',fontSize,'FontWeight','bold'); ylabel('[m]','FontSize',fontSize,'FontWeight','bold'); CH = colorbar; CH.Label.String = 'W/m^{2}'; axis square; set(gca,'FontSize',fontSize,'FontWeight','normal');
axes(ha(5)); imagesc(xCoords(dispRngA{1}),yCoords(dispRngA{2}),compIntensity(PSFlensPlaneDistribution(dispRngA{1},dispRngA{2}),irradianceScaling)); title(sprintf('Plane at L3: EFL=%1.1fm',focalLengthB)); xlabel('[m]','FontSize',fontSize,'FontWeight','bold'); ylabel('[m]','FontSize',fontSize,'FontWeight','bold'); CH = colorbar; CH.Label.String = 'W/m^{2}'; axis square; set(gca,'FontSize',fontSize,'FontWeight','normal');
axes(ha(6)); imagesc(xCoords(dispRngB{1}),yCoords(dispRngB{2}),outputIntensity(dispRngB{1},dispRngB{2})); title(sprintf('PSF Distribution: z=%1.1fm',3*focalLengthA + 2*focalLengthB)); xlabel('[m]','FontSize',fontSize,'FontWeight','bold'); ylabel('[m]','FontSize',fontSize,'FontWeight','bold'); CH = colorbar; CH.Label.String = 'W/m^{2}'; axis square; set(gca,'FontSize',fontSize,'FontWeight','normal');
set(figureHandleB,'units','normalized','position',[0,0,1,1]);

figureHandleC = figure('Color', 'white');
[ha, ~] = tight_subplot(2,3,plotGaps,heighMargins,widthMargins);
axes(ha(1)); imagesc(xCoords(dispRngA{1}),yCoords(dispRngA{2}),compAngle(inputPlane(dispRngA{1},dispRngA{2}))); title('Input Plane','FontSize',fontSize,'FontWeight','bold'); xlabel('[m]','FontSize',fontSize,'FontWeight','bold'); ylabel('[m]','FontSize',fontSize,'FontWeight','bold'); CH = colorbar; CH.Label.String = 'rad'; axis square; set(gca,'FontSize',fontSize,'FontWeight','normal');
axes(ha(2)); imagesc(xCoords(dispRngA{1}),yCoords(dispRngA{2}),compAngle(focalPlaneDistribution(dispRngA{1},dispRngA{2}))); title(sprintf('Focal Plane: EFL=%1.1fm',focalLengthA)); xlabel('[m]','FontSize',fontSize,'FontWeight','bold'); ylabel('[m]','FontSize',fontSize,'FontWeight','bold'); CH = colorbar; CH.Label.String = 'rad'; axis square; set(gca,'FontSize',fontSize,'FontWeight','normal');
axes(ha(3)); imagesc(xCoords(dispRngA{1}),yCoords(dispRngA{2}),compAngle(ocularPlaneDistribution(dispRngA{1},dispRngA{2}))); title(sprintf('Plane at L2: z=%1.1fm',2*focalLengthA)); xlabel('[m]','FontSize',fontSize,'FontWeight','bold'); ylabel('[m]','FontSize',fontSize,'FontWeight','bold'); CH = colorbar; CH.Label.String = 'rad'; axis square; set(gca,'FontSize',fontSize,'FontWeight','normal');
axes(ha(4)); imagesc(xCoords(dispRngA{1}),yCoords(dispRngA{2}),compAngle(truncatedLyotPlane(dispRngA{1},dispRngA{2}))); title(sprintf('Truncated Lyot: %1.1fmm', LyotApertureDiameter*10^abs(ceil(log10(LyotApertureDiameter)-1)))); xlabel('[m]','FontSize',fontSize,'FontWeight','bold'); ylabel('[m]','FontSize',fontSize,'FontWeight','bold'); CH = colorbar; CH.Label.String = 'rad'; axis square; set(gca,'FontSize',fontSize,'FontWeight','normal');
axes(ha(5)); imagesc(xCoords(dispRngA{1}),yCoords(dispRngA{2}),compAngle(PSFlensPlaneDistribution(dispRngA{1},dispRngA{2}))); title(sprintf('Plane at L3: EFL=%1.1fm',focalLengthB)); xlabel('[m]','FontSize',fontSize,'FontWeight','bold'); ylabel('[m]','FontSize',fontSize,'FontWeight','bold'); CH = colorbar; CH.Label.String = 'rad'; axis square; set(gca,'FontSize',fontSize,'FontWeight','normal');
axes(ha(6)); imagesc(xCoords(dispRngA{1}),yCoords(dispRngA{2}),compAngle(PSFplaneDistribution(dispRngA{1},dispRngA{2}))); title(sprintf('PSF Distribution: z=%1.1fm',3*focalLengthA + 2*focalLengthB)); xlabel('[m]','FontSize',fontSize,'FontWeight','bold'); ylabel('[m]','FontSize',fontSize,'FontWeight','bold'); CH = colorbar; CH.Label.String = 'rad'; axis square; set(gca,'FontSize',fontSize,'FontWeight','normal');
set(figureHandleC,'units','normalized','position',[0,0,1,1]);