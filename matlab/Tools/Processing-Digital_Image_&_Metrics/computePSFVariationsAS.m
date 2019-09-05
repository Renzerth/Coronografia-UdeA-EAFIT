function [LyotPlaneIntensities, PSFplaneIntensities, LyotReference, PSFreference] = computePSFVariationsAS(TC, glvect)
%% Field Properties definition
compIntensity = @(complexField, energyScaling) energyScaling*abs(complexField).^2;

%% Plane properties
planeSize = 2*[2.54e-2, 2.54e-2]; % m
spaceSamples = [1024, 1024];
halfSize = floor((spaceSamples+1)/2);
shiftDistance = 0;

%% Filter Properties
tcvect = 1:TC;
grayLevels = length(glvect);
totalFields = TC*grayLevels;

%% Data Initialization
LyotPlaneIntensities = cell(1,totalFields);
PSFplaneIntensities = cell(1,totalFields);

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

%% Spatial coordinates
xCenter= floor((spaceSamples(1) + 1)/2);
yCenter = floor((spaceSamples(2) + 1)/2);
[~, ~, ~, ~, theta, rho] = makeSpaceCoordinates(planeSize(1),planeSize(2),spaceSamples(1),spaceSamples(2),xCenter,yCenter);

%% Spectral coordinates
samplingFactor = 1;
[freqVectX, freqVectY,~,~,analysisScaling,normNMFactor,synthesisScaling] = computeFreqVector(planeSize, spaceSamples, samplingFactor);

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
zernikeCoeffs = zeros(1,10);
zernikeCoeffs(aberrationIndex) = aberrationValue;
[systemPhase,~] = Zernike_Builder(zernikeCoeffs',aberrationPupilRadius/planeSize(1),min(spaceSamples),false);
systemPhase(isnan(systemPhase)) = 0;

%% Optical System Phase Properties Definition
Pupil = double(rho <= aberrationPupilRadius);
pupilTransferFunct = Pupil.*exp(1i*systemPhase);

%% Compute Propagators
propDistance = focalLengthA - shiftDistance;
[dataCoordX, dataCoordY, freqMeshX, freqMeshY, ~, freqCirc] = prepareSpectralProp_mod(spaceSamples, freqVectX, freqVectY, lambda);
distanceShiftA = exp(-1i*k*focalLengthA);
distanceShiftB = exp(-1i*k*focalLengthB);

[focalPlanePropagationKernelA] = freqCirc.*exp(1i*k*focalLengthA*sqrt(1 - (lambda*freqMeshX).^2 - (lambda*freqMeshY).^2));
[focalPlanePropagationKernelB] = freqCirc.*exp(1i*k*focalLengthB*sqrt(1 - (lambda*freqMeshX).^2 - (lambda*freqMeshY).^2));
[shiftDistanceKernel] = freqCirc.*exp(1i*k*propDistance*sqrt(1 - (lambda*freqMeshX).^2 - (lambda*freqMeshY).^2));

%% Generate Focal Distribution
systemTransFunct = pupilTransferFunct.*shiftDistanceKernel;
[focalPlaneDistribution] = distanceShiftA*convoluteSignal(systemTransFunct, inputPlane, normNMFactor, analysisScaling, synthesisScaling, dataCoordX, dataCoordY, halfSize);

%% Lyot Aperture Definition
LyotAperture = double(rho < LyotApertureDiameter/2);

%% Coronograph End-To-End Propagation with Dynamic Filter
for tcIndex = 1:TC
    for glIndex = 1:grayLevels
        
        singleIndex = glIndex + (tcIndex - 1)*grayLevels;
        fprintf('Processing Iteration %d/%d.\n\r', singleIndex, totalFields);
        
%         vortexMask = exp(1i*(tcvect(tcIndex)*(theta - pi/2)));
        vortexMask = spiralGen2(spaceSamples,tcvect(tcIndex));
        [discretMap] = discretizeMap(angle(vortexMask),glvect(glIndex));
        vortexMask = exp(1i*discretMap);
        
        SLMfilteredField = focalPlaneDistribution.*vortexMask;
        [ocularPlaneDistribution] = distanceShiftA*deconvoluteSignal(focalPlanePropagationKernelA, SLMfilteredField, normNMFactor, analysisScaling, synthesisScaling, dataCoordX, dataCoordY, halfSize);
        LyotPlaneDistribution = ocularPlaneDistribution.*diverLens;
        
        truncatedLyotPlane = distanceShiftA*LyotPlaneDistribution.*LyotAperture;
        PSFlensPlaneDistribution = distanceShiftB*truncatedLyotPlane.*lensB;
        [PSFplaneDistribution] = distanceShiftB*convoluteSignal(focalPlanePropagationKernelB, PSFlensPlaneDistribution, normNMFactor, analysisScaling, synthesisScaling, dataCoordX, dataCoordY, halfSize);
        
        LyotPlaneIntensities{1, singleIndex} = compIntensity(LyotPlaneDistribution, irradianceScaling);
        PSFplaneIntensities{1, singleIndex} = compIntensity(PSFplaneDistribution, irradianceScaling);
        
    end
end
%% Compute PSF and Lyot References
[ocularPlaneDistribution] = distanceShiftA*deconvoluteSignal(focalPlanePropagationKernelA, focalPlaneDistribution, normNMFactor, analysisScaling, synthesisScaling, dataCoordX, dataCoordY, halfSize);
LyotPlaneDistribution = ocularPlaneDistribution.*diverLens;

truncatedLyotPlane = distanceShiftA*LyotPlaneDistribution.*LyotAperture;
PSFlensPlaneDistribution = distanceShiftB*truncatedLyotPlane.*lensB;
[PSFplaneDistribution] = distanceShiftB*convoluteSignal(focalPlanePropagationKernelB, PSFlensPlaneDistribution, normNMFactor, analysisScaling, synthesisScaling, dataCoordX, dataCoordY, halfSize);

PSFreference = compIntensity(PSFplaneDistribution, irradianceScaling);
LyotReference = compIntensity(LyotPlaneDistribution, irradianceScaling);
end