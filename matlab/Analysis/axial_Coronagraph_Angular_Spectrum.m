%% path dependencies
addpath(genpath(fileparts(pwd)));

%% Field Properties definition
fieldIntensity = @(complexField, energyScaling) energyScaling*abs(complexField).^2;
fieldAngle = @(complexField) angle(complexField);

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
irradianceScaling = insidentEnergy*(propagationSpeed*mediumRefracIndex*mediumElecPermitivity)/2; % [W/m^2]

%% System Properties
illuminationDiameter = 3e-3; %[m]
obstructionDiameter = illuminationDiameter/2;
LyotApertureDiameter = illuminationDiameter; % [m]
aberrationPupilRadius = max(planeSize); %[m]

focalLengthA = 2; % [m]
focalLengthB = 2;
aberrationIndex = 8;
aberrationValue = 0*pi*(1.0021)*0.01;
FresnelNumber = (illuminationDiameter)^2/(focalLengthA*lambda);

%% Spatial coordinates
xCenter = floor((spaceSamples(1) + 1)/2);
yCenter = floor((spaceSamples(2) + 1)/2);
[xCoords, yCoords, X, Y, theta, rho] = makeSpaceCoordinates(planeSize(1),planeSize(2),spaceSamples(1),spaceSamples(2),xCenter,yCenter);

%% Spectral coordinates
samplingFactor = 1;
[freqVectX,freqVectY,spXperiod,spYperiod,analysisScaling,normNMFactor,synthesisScaling] = computeFreqVector(planeSize, spaceSamples, samplingFactor);

%% Vortex Mask properties
TC = 2;
grayLevels = 24;
vortexMask = spiralGen2(spaceSamples,TC);
[vortexMask] = discretizeMap(angle(vortexMask),grayLevels);

%% Lens properties
lensRadii = rho;
lensDiameter = 2.54e-2;
lensAperture = double(rho <= lensDiameter/2);
lensA = lensAperture.*exp(-1i*k/(2*focalLengthA)*(lensRadii).^2);
lensB = lensAperture.*exp(-1i*k/(2*focalLengthB)*(lensRadii).^2);
diverLens = lensAperture.*exp(1i*k/(2*focalLengthA)*(lensRadii).^2);

%% Uniform light definition after first focal lens
% aperture = double(rho <= illuminationDiameter/2).*double(rho >= obstructionDiameter/2);
aperture = double(rho <= illuminationDiameter/2);
inputPlane =  insidentEnergy*aperture;
% inputPlane = insidentEnergy*aperture.*evaluateGaussianField(X,Y,illuminationDiameter/4,[0,0],lambda);

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

angularSpectrum = sqrt(1 - (lambda*freqMeshX).^2 - (lambda*freqMeshY).^2);
propagationKernel = @ (freqLowPass, k, angularSpectrum, z) freqLowPass.*exp(1i*k*z*angularSpectrum);

midPlanesNumber = 20;
propagationDistancesA = linspace(1e-5,focalLengthA,midPlanesNumber);
propagationDistancesB = linspace(1e-5,focalLengthB,midPlanesNumber);
extendedAxialDistance = linspace(0,6*focalLengthB,6*midPlanesNumber);

%% Axial analysis
windowFraction = 5;
windowSize = floor(spaceSamples(2)/windowFraction);
halfwindow = floor(windowSize/2);
transverseCutRange =  halfSize(2) - halfwindow:halfSize(2) + halfwindow;  % centered window of quarter of size
cutDataLength = length(transverseCutRange);

interpFactor = 2;
interpRangeX = 1:1/interpFactor^2:midPlanesNumber;
interpRangeY = 1:1/interpFactor^2:cutDataLength;

interpDIstanceA = interp1(propagationDistancesA,interpRangeX);
interpTransverseCoords = interp1(yCoords(transverseCutRange),interpRangeY);

interpRangeExtendedX = 1:1/interpFactor^2:6*midPlanesNumber;
interpExtendedDIstanceA = interp1(extendedAxialDistance,interpRangeExtendedX);

%% Entrance Propagation
targetPlane = inputPlane;
[entrancePropagation,entranceCut] = propagatePerPlane(targetPlane, k, angularSpectrum, propagationKernel, spaceSamples, halfSize,transverseCutRange,cutDataLength, propagationDistancesA,midPlanesNumber,normNMFactor, analysisScaling, synthesisScaling, dataCoordX, dataCoordY, freqCirc);

%% L1 focalizing Propagation
targetPlane = entrancePropagation(:,:,end).*lensA;
[focalizingL1Propagation,focalizingL1Cut] = propagatePerPlane(targetPlane, k, angularSpectrum, propagationKernel, spaceSamples, halfSize,transverseCutRange,cutDataLength, propagationDistancesA,midPlanesNumber,normNMFactor, analysisScaling, synthesisScaling, dataCoordX, dataCoordY, freqCirc);

%% Diverging Propagation
targetPlane = focalizingL1Propagation(:,:,end).*vortexMask;
[divergingL1Propagation,divergingL1Cut] = propagatePerPlane(targetPlane, k, angularSpectrum, propagationKernel, spaceSamples, halfSize,transverseCutRange,cutDataLength, propagationDistancesA+0.01,midPlanesNumber,normNMFactor, analysisScaling, synthesisScaling, dataCoordX, dataCoordY, freqCirc);

%% Generate Lyot distribution - Pupil Plane
targetPlane = divergingL1Propagation(:,:,end).*lensA;
[planeL2Propagation,planeL2Cut] = propagatePerPlane(targetPlane, k, angularSpectrum, propagationKernel, spaceSamples, halfSize,transverseCutRange,cutDataLength, propagationDistancesA,midPlanesNumber,normNMFactor, analysisScaling, synthesisScaling, dataCoordX, dataCoordY, freqCirc);

%% Lyot Pupil Truncation Propagation
LyotAperture = double(rho < LyotApertureDiameter/2);
targetPlane = planeL2Propagation(:,:,end).*LyotAperture;
[planeLyotPropagation,planeLyotCut] = propagatePerPlane(targetPlane, k, angularSpectrum, propagationKernel, spaceSamples, halfSize,transverseCutRange,cutDataLength, propagationDistancesA,midPlanesNumber,normNMFactor, analysisScaling, synthesisScaling, dataCoordX, dataCoordY, freqCirc);

%% PSF Propagation
targetPlane = planeLyotPropagation(:,:,end).*lensB;
[focalizingL3Propagation,focalizingL3Cut] = propagatePerPlane(targetPlane, k, angularSpectrum, propagationKernel, spaceSamples, halfSize,transverseCutRange,cutDataLength, propagationDistancesB,midPlanesNumber,normNMFactor, analysisScaling, synthesisScaling, dataCoordX, dataCoordY, freqCirc);

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

%% Visualize Coronagraphic Sequence -- Intensities
makeGifEnabled = false
if makeGifEnabled
    timeDelay = 0.2;
    animationDelay = 0.3;
    intensityScale = 1;

    gifNameIntensities = 'axial_coronagraph_intensities.gif';
    gifNamePhases = 'axial_coronagraph_phases.gif';

    figure('Color', 'white');
    updateDisplayHandler = imagesc(xCoords(dispRngA{1}),yCoords(dispRngA{2}),fieldIntensity(entrancePropagation(dispRngA{1},dispRngA{2},1),intensityScale)); 
    xlabel('[m]','FontSize',fontSize,'FontWeight','bold');
    ylabel('[m]','FontSize',fontSize,'FontWeight','bold');
    set(gca,'FontSize',fontSize,'FontWeight','normal');
    grid on; colorbar; grid on; axis square; set(gca,'GridColor',[1,1,1]);

    for imageIndex = 1:midPlanesNumber ; title(sprintf('z=%1.1fm',propagationDistancesA(imageIndex))); set(updateDisplayHandler,'CData',fieldIntensity(entrancePropagation(dispRngA{1},dispRngA{2},imageIndex),intensityScale)); colorbar; axis square; Make_Gif(gifNameIntensities,imageIndex,animationDelay); pause(timeDelay); end
    for imageIndex = 1:midPlanesNumber ; title(sprintf('z=%1.1fm',1*focalLengthA + propagationDistancesA(imageIndex)));set(updateDisplayHandler,'CData',fieldIntensity(focalizingL1Propagation(dispRngA{1},dispRngA{2},imageIndex),intensityScale)); Make_Gif(gifNameIntensities,midPlanesNumber+imageIndex,animationDelay); colorbar; axis square; pause(timeDelay); end
    for imageIndex = 1:midPlanesNumber ; title(sprintf('z=%1.1fm',2*focalLengthA + propagationDistancesA(imageIndex)));set(updateDisplayHandler,'CData',fieldIntensity(divergingL1Propagation(dispRngA{1},dispRngA{2},imageIndex),intensityScale)); Make_Gif(gifNameIntensities,2*midPlanesNumber+imageIndex,animationDelay); colorbar; axis square; pause(timeDelay); end
    for imageIndex = 1:midPlanesNumber ; title(sprintf('z=%1.1fm',3*focalLengthA + propagationDistancesA(imageIndex)));set(updateDisplayHandler,'CData',fieldIntensity(planeL2Propagation(dispRngA{1},dispRngA{2},imageIndex),intensityScale)); Make_Gif(gifNameIntensities,3*midPlanesNumber+imageIndex,animationDelay); colorbar; axis square;  pause(timeDelay); end
    for imageIndex = 1:midPlanesNumber ; title(sprintf('z=%1.1fm',4*focalLengthA + propagationDistancesA(imageIndex)));set(updateDisplayHandler,'CData',fieldIntensity(planeLyotPropagation(dispRngA{1},dispRngA{2},imageIndex),intensityScale)); Make_Gif(gifNameIntensities,4*midPlanesNumber+imageIndex,animationDelay); colorbar; axis square; pause(timeDelay); end
    for imageIndex = 1:midPlanesNumber ; title(sprintf('z=%1.1fm',5*focalLengthB + propagationDistancesB(imageIndex)));set(updateDisplayHandler,'CData',fieldIntensity(focalizingL3Propagation(dispRngA{1},dispRngA{2},imageIndex),intensityScale)); Make_Gif(gifNameIntensities,5*midPlanesNumber+imageIndex,animationDelay); colorbar; axis square; pause(timeDelay); end

    %% Visualize Coronagraphic Sequence -- Phases
    close(gcf);
    figure('Color', 'white');
    updateDisplayHandler = imagesc(xCoords(dispRngA{1}),yCoords(dispRngA{2}),fieldAngle(entrancePropagation(dispRngA{1},dispRngA{2},1))); 
    xlabel('[m]','FontSize',fontSize,'FontWeight','bold');
    ylabel('[m]','FontSize',fontSize,'FontWeight','bold');
    set(gca,'FontSize',fontSize,'FontWeight','normal');
    grid on; colorbar; grid on; axis square; set(gca,'GridColor',[1,1,1]);

    for imageIndex = 1:midPlanesNumber ; title(sprintf('z=%1.1fm',propagationDistancesA(imageIndex))); set(updateDisplayHandler,'CData',fieldAngle(entrancePropagation(dispRngA{1},dispRngA{2},imageIndex))); colorbar; axis square; Make_Gif(gifNamePhases,imageIndex,animationDelay); pause(timeDelay); end
    for imageIndex = 1:midPlanesNumber ; title(sprintf('z=%1.1fm',1*focalLengthA + propagationDistancesA(imageIndex)));set(updateDisplayHandler,'CData',fieldAngle(focalizingL1Propagation(dispRngA{1},dispRngA{2},imageIndex))); Make_Gif(gifNamePhases,midPlanesNumber+imageIndex,animationDelay); colorbar; axis square; pause(timeDelay); end
    for imageIndex = 1:midPlanesNumber ; title(sprintf('z=%1.1fm',2*focalLengthA + propagationDistancesA(imageIndex)));set(updateDisplayHandler,'CData',fieldAngle(divergingL1Propagation(dispRngA{1},dispRngA{2},imageIndex))); Make_Gif(gifNamePhases,2*midPlanesNumber+imageIndex,animationDelay); colorbar; axis square; pause(timeDelay); end
    for imageIndex = 1:midPlanesNumber ; title(sprintf('z=%1.1fm',3*focalLengthA + propagationDistancesA(imageIndex)));set(updateDisplayHandler,'CData',fieldAngle(planeL2Propagation(dispRngA{1},dispRngA{2},imageIndex))); Make_Gif(gifNamePhases,3*midPlanesNumber+imageIndex,animationDelay); colorbar; axis square;  pause(timeDelay); end
    for imageIndex = 1:midPlanesNumber ; title(sprintf('z=%1.1fm',4*focalLengthA + propagationDistancesA(imageIndex)));set(updateDisplayHandler,'CData',fieldAngle(planeLyotPropagation(dispRngA{1},dispRngA{2},imageIndex))); Make_Gif(gifNamePhases,4*midPlanesNumber+imageIndex,animationDelay); colorbar; axis square; pause(timeDelay); end
    for imageIndex = 1:midPlanesNumber ; title(sprintf('z=%1.1fm',5*focalLengthB + propagationDistancesB(imageIndex)));set(updateDisplayHandler,'CData',fieldAngle(focalizingL3Propagation(dispRngA{1},dispRngA{2},imageIndex))); Make_Gif(gifNamePhases,5*midPlanesNumber+imageIndex,animationDelay); colorbar; axis square; pause(timeDelay); end
end
%% Axial Distribution Behavior
figure('color','white'); imagesc(interpDIstanceA,interpTransverseCoords,interp2(entranceCut,interpFactor))
title('Input Plane Axial Propagation - 1F','FontSize',fontSize,'FontWeight','bold')
xlabel('Axial Distance [m]','FontSize',fontSize,'FontWeight','bold');
ylabel('Transverse Cut [m]','FontSize',fontSize,'FontWeight','bold');
set(gca,'FontSize',fontSize,'FontWeight','normal'); grid on
colorbar; grid on; axis tight; set(gca,'GridColor',[1,1,1]);

figure('color','white'); imagesc(interpDIstanceA,interpTransverseCoords,interp2(focalizingL1Cut,interpFactor))
title('L1 Axial Focalization - 2F','FontSize',fontSize,'FontWeight','bold')
xlabel('Axial Distance [m]','FontSize',fontSize,'FontWeight','bold');
ylabel('Transverse Cut [m]','FontSize',fontSize,'FontWeight','bold');
set(gca,'FontSize',fontSize,'FontWeight','normal'); grid on
colorbar; grid on; axis tight; set(gca,'GridColor',[1,1,1]);

figure('color','white'); imagesc(interpDIstanceA,interpTransverseCoords,interp2(divergingL1Cut,interpFactor))
title('Diverging Plane to L2 - 3F','FontSize',fontSize,'FontWeight','bold')
xlabel('Axial Distance [m]','FontSize',fontSize,'FontWeight','bold');
ylabel('Transverse Cut [m]','FontSize',fontSize,'FontWeight','bold');
set(gca,'FontSize',fontSize,'FontWeight','normal'); grid on
colorbar; grid on; axis tight; set(gca,'GridColor',[1,1,1]);

figure('color','white'); imagesc(interpDIstanceA,interpTransverseCoords,interp2(planeL2Cut,interpFactor))
title('L2 Axial Focalization - 4F','FontSize',fontSize,'FontWeight','bold')
xlabel('Axial Distance [m]','FontSize',fontSize,'FontWeight','bold');
ylabel('Transverse Cut [m]','FontSize',fontSize,'FontWeight','bold');
set(gca,'FontSize',fontSize,'FontWeight','normal'); grid on
colorbar; grid on; axis tight; set(gca,'GridColor',[1,1,1]);

figure('color','white'); imagesc(interpDIstanceA,interpTransverseCoords,interp2(planeLyotCut,interpFactor))
title('Diverging Plane to L3 - 5F','FontSize',fontSize,'FontWeight','bold')
xlabel('Axial Distance [m]','FontSize',fontSize,'FontWeight','bold');
ylabel('Transverse Cut [m]','FontSize',fontSize,'FontWeight','bold');
set(gca,'FontSize',fontSize,'FontWeight','normal'); grid on
colorbar; grid on; axis tight; set(gca,'GridColor',[1,1,1]);

figure('color','white'); imagesc(interpDIstanceA,interpTransverseCoords,interp2(focalizingL3Cut,interpFactor))
title('L3 Axial Focalization - 6F','FontSize',fontSize,'FontWeight','bold')
xlabel('Axial Distance [m]','FontSize',fontSize,'FontWeight','bold');
ylabel('Transverse Cut [m]','FontSize',fontSize,'FontWeight','bold');
set(gca,'FontSize',fontSize,'FontWeight','normal'); grid on
colorbar; grid on; axis tight; set(gca,'GridColor',[1,1,1]);

%% Comple Axial Plots
completeSetNorm = interp2([bsxfun(@rdivide,entranceCut,max(entranceCut)),bsxfun(@rdivide,focalizingL1Cut,max(focalizingL1Cut)),bsxfun(@rdivide,divergingL1Cut,max(divergingL1Cut)),bsxfun(@rdivide,planeL2Cut,max(planeL2Cut)),bsxfun(@rdivide,planeLyotCut,max(planeLyotCut)),bsxfun(@rdivide,focalizingL3Cut,max(focalizingL3Cut))], interpFactor);
figure('color','white'); imagesc(interpExtendedDIstanceA, interpTransverseCoords, completeSetNorm);
title('Complete Axial propagation 1F - 6F: Normalized','FontSize',fontSize,'FontWeight','bold')
xlabel('Axial Distance [m]','FontSize',fontSize,'FontWeight','bold');
ylabel('Transverse Cut [m]','FontSize',fontSize,'FontWeight','bold');
set(gca,'FontSize',fontSize,'FontWeight','normal'); grid on
colorbar; grid on; axis tight; set(gca,'GridColor',[1,1,1]);

completeSet = interp2([entranceCut,focalizingL1Cut,divergingL1Cut,planeL2Cut,planeLyotCut,focalizingL3Cut], interpFactor);
figure('color','white'); imagesc(interpExtendedDIstanceA, interpTransverseCoords,10*log10(completeSet));
title('Complete Axial propagation: LogView','FontSize',fontSize,'FontWeight','bold')
xlabel('Axial Distance [m]','FontSize',fontSize,'FontWeight','bold');
ylabel('Transverse Cut [m]','FontSize',fontSize,'FontWeight','bold');
set(gca,'FontSize',fontSize,'FontWeight','normal'); grid on
grid on; axis tight; set(gca,'GridColor',[0,0,0]); colormap('gray')
CH = colorbar; CH.Label.String = 'dB'; set(gca,'FontSize',fontSize,'FontWeight','normal');
set(gca,'xticklabels',{'0F','1F','2F','3F','4F','5F','6F'});

%% Relevant Planes
focalPlanePropagationKernelA = propagationKernel(freqCirc, k, angularSpectrum, propagationDistancesA(midPlanesNumber));

inputIntensity = fieldIntensity(entrancePropagation(:,:,1),irradianceScaling);
focalPlaneDistribution = focalizingL1Propagation(:,:,end);
ocularPlaneDistribution = divergingL1Propagation(:,:,end).*lensA;
truncatedLyotPlane = planeL2Propagation(:,:,end).*LyotAperture;
PSFlensPlaneDistribution = planeLyotPropagation(:,:,end).*lensB;
PSFplaneDistribution = focalizingL3Propagation(:,:,end);
outputIntensity = fieldIntensity(PSFplaneDistribution,irradianceScaling);

%% Plots

figureHandleA = figure('Color', 'white');
[ha, ~] = tight_subplot(2,3,plotGaps,heighMargins,widthMargins);
axes(ha(1)); imagesc(xCoords,yCoords,inputIntensity); title('Input Plane'); xlabel('[m]','FontSize',fontSize,'FontWeight','bold'); ylabel('[m]','FontSize',fontSize,'FontWeight','bold'); CH = colorbar; CH.Label.String = 'W/m^{2}'; axis square; set(gca,'FontSize',fontSize,'FontWeight','normal');
axes(ha(2)); imagesc(xCoords,yCoords,outputIntensity); title('Output Plane'); xlabel('[m]','FontSize',fontSize,'FontWeight','bold'); ylabel('[m]','FontSize',fontSize,'FontWeight','bold'); CH = colorbar; CH.Label.String = 'W/m^{2}'; axis square; set(gca,'FontSize',fontSize,'FontWeight','normal');
axes(ha(3)); imagesc(xCoords,yCoords,angle(lensA)); title(sprintf('L1 Lens Phase: f=%1.1fm',focalLengthA)); xlabel('[m]','FontSize',fontSize,'FontWeight','bold'); ylabel('[m]','FontSize',fontSize,'FontWeight','bold'); CH = colorbar; CH.Label.String = 'rad'; axis square; set(gca,'FontSize',fontSize,'FontWeight','normal');
axes(ha(4)); imagesc(xCoords,yCoords,10*log10(inputIntensity)); title('Input Log10 view'); xlabel('[m]','FontSize',fontSize,'FontWeight','bold'); ylabel('[m]','FontSize',fontSize,'FontWeight','bold'); CH = colorbar; CH.Label.String = 'dB'; axis square; set(gca,'FontSize',fontSize,'FontWeight','normal');
axes(ha(5)); imagesc(xCoords,yCoords,10*log10(outputIntensity)); title('Output Log10 view'); xlabel('[m]','FontSize',fontSize,'FontWeight','bold'); ylabel('[m]','FontSize',fontSize,'FontWeight','bold');CH = colorbar; CH.Label.String = 'dB'; axis square; set(gca,'FontSize',fontSize,'FontWeight','normal');
axes(ha(6)); imagesc(freqVectX,freqVectY,angle(focalPlanePropagationKernelA)); title(sprintf('Vaccum Transfer Function: D=%1.1fm', propDistance)); xlabel('fx [1/m]','FontSize',fontSize,'FontWeight','bold'); ylabel('fy [1/m]','FontSize',fontSize,'FontWeight','bold'); CH = colorbar; CH.Label.String = 'rad'; axis square; set(gca,'FontSize',fontSize,'FontWeight','normal');
set(figureHandleA,'units','normalized','position',[0,0,1,1]);

figureHandleB = figure('Color', 'white');
[ha, ~] = tight_subplot(2,3,plotGaps,heighMargins,widthMargins);
axes(ha(1)); imagesc(xCoords(dispRngA{1}),yCoords(dispRngA{2}),inputIntensity(dispRngA{1},dispRngA{2})); title('Input Plane','FontSize',fontSize,'FontWeight','bold'); xlabel('[m]','FontSize',fontSize,'FontWeight','bold'); ylabel('[m]','FontSize',fontSize,'FontWeight','bold'); CH = colorbar; CH.Label.String = 'W/m^{2}'; axis square; set(gca,'FontSize',fontSize,'FontWeight','normal');
axes(ha(2)); imagesc(xCoords(dispRngA{1}),yCoords(dispRngA{2}),fieldIntensity(focalPlaneDistribution(dispRngA{1},dispRngA{2}),irradianceScaling)); title(sprintf('Focal Plane: EFL=%1.1fm',focalLengthA)); xlabel('[m]','FontSize',fontSize,'FontWeight','bold'); ylabel('[m]','FontSize',fontSize,'FontWeight','bold'); CH = colorbar; CH.Label.String = 'W/m^{2}'; axis square; set(gca,'FontSize',fontSize,'FontWeight','normal');
axes(ha(3)); imagesc(xCoords(dispRngA{1}),yCoords(dispRngA{2}),fieldIntensity(ocularPlaneDistribution(dispRngA{1},dispRngA{2}),irradianceScaling)); title(sprintf('Plane at L2: z=%1.1fm',2*focalLengthA)); xlabel('[m]','FontSize',fontSize,'FontWeight','bold'); ylabel('[m]','FontSize',fontSize,'FontWeight','bold'); CH = colorbar; CH.Label.String = 'W/m^{2}'; axis square; set(gca,'FontSize',fontSize,'FontWeight','normal');
axes(ha(4)); imagesc(xCoords(dispRngA{1}),yCoords(dispRngA{2}),fieldIntensity(truncatedLyotPlane(dispRngA{1},dispRngA{2}),irradianceScaling)); title(sprintf('Truncated Lyot: %1.1fmm', LyotApertureDiameter*10^abs(ceil(log10(LyotApertureDiameter)-1)))); xlabel('[m]','FontSize',fontSize,'FontWeight','bold'); ylabel('[m]','FontSize',fontSize,'FontWeight','bold'); CH = colorbar; CH.Label.String = 'W/m^{2}'; axis square; set(gca,'FontSize',fontSize,'FontWeight','normal');
axes(ha(5)); imagesc(xCoords(dispRngA{1}),yCoords(dispRngA{2}),fieldIntensity(PSFlensPlaneDistribution(dispRngA{1},dispRngA{2}),irradianceScaling)); title(sprintf('Plane at L3: EFL=%1.1fm',focalLengthB)); xlabel('[m]','FontSize',fontSize,'FontWeight','bold'); ylabel('[m]','FontSize',fontSize,'FontWeight','bold'); CH = colorbar; CH.Label.String = 'W/m^{2}'; axis square; set(gca,'FontSize',fontSize,'FontWeight','normal');
axes(ha(6)); imagesc(xCoords(dispRngB{1}),yCoords(dispRngB{2}),outputIntensity(dispRngB{1},dispRngB{2})); title(sprintf('PSF Distribution: z=%1.1fm',4*focalLengthA + 2*focalLengthB)); xlabel('[m]','FontSize',fontSize,'FontWeight','bold'); ylabel('[m]','FontSize',fontSize,'FontWeight','bold'); CH = colorbar; CH.Label.String = 'W/m^{2}'; axis square; set(gca,'FontSize',fontSize,'FontWeight','normal');
set(figureHandleB,'units','normalized','position',[0,0,1,1]);

figureHandleC = figure('Color', 'white');
[ha, ~] = tight_subplot(2,3,plotGaps,heighMargins,widthMargins);
axes(ha(1)); imagesc(xCoords(dispRngA{1}),yCoords(dispRngA{2}),fieldAngle(inputPlane(dispRngA{1},dispRngA{2}))); title('Input Plane','FontSize',fontSize,'FontWeight','bold'); xlabel('[m]','FontSize',fontSize,'FontWeight','bold'); ylabel('[m]','FontSize',fontSize,'FontWeight','bold'); CH = colorbar; CH.Label.String = 'rad'; axis square; set(gca,'FontSize',fontSize,'FontWeight','normal');
axes(ha(2)); imagesc(xCoords(dispRngA{1}),yCoords(dispRngA{2}),fieldAngle(focalPlaneDistribution(dispRngA{1},dispRngA{2}))); title(sprintf('Focal Plane: EFL=%1.1fm',focalLengthA)); xlabel('[m]','FontSize',fontSize,'FontWeight','bold'); ylabel('[m]','FontSize',fontSize,'FontWeight','bold'); CH = colorbar; CH.Label.String = 'rad'; axis square; set(gca,'FontSize',fontSize,'FontWeight','normal');
axes(ha(3)); imagesc(xCoords(dispRngA{1}),yCoords(dispRngA{2}),fieldAngle(ocularPlaneDistribution(dispRngA{1},dispRngA{2}))); title(sprintf('Plane at L2: z=%1.1fm',2*focalLengthA)); xlabel('[m]','FontSize',fontSize,'FontWeight','bold'); ylabel('[m]','FontSize',fontSize,'FontWeight','bold'); CH = colorbar; CH.Label.String = 'rad'; axis square; set(gca,'FontSize',fontSize,'FontWeight','normal');
axes(ha(4)); imagesc(xCoords(dispRngA{1}),yCoords(dispRngA{2}),fieldAngle(truncatedLyotPlane(dispRngA{1},dispRngA{2}))); title(sprintf('Truncated Lyot: %1.1fmm', LyotApertureDiameter*10^abs(ceil(log10(LyotApertureDiameter)-1)))); xlabel('[m]','FontSize',fontSize,'FontWeight','bold'); ylabel('[m]','FontSize',fontSize,'FontWeight','bold'); CH = colorbar; CH.Label.String = 'rad'; axis square; set(gca,'FontSize',fontSize,'FontWeight','normal');
axes(ha(5)); imagesc(xCoords(dispRngA{1}),yCoords(dispRngA{2}),fieldAngle(PSFlensPlaneDistribution(dispRngA{1},dispRngA{2}))); title(sprintf('Plane at L3: EFL=%1.1fm',focalLengthB)); xlabel('[m]','FontSize',fontSize,'FontWeight','bold'); ylabel('[m]','FontSize',fontSize,'FontWeight','bold'); CH = colorbar; CH.Label.String = 'rad'; axis square; set(gca,'FontSize',fontSize,'FontWeight','normal');
axes(ha(6)); imagesc(xCoords(dispRngA{1}),yCoords(dispRngA{2}),fieldAngle(PSFplaneDistribution(dispRngA{1},dispRngA{2}))); title(sprintf('PSF Distribution: z=%1.1fm',4*focalLengthA + 2*focalLengthB)); xlabel('[m]','FontSize',fontSize,'FontWeight','bold'); ylabel('[m]','FontSize',fontSize,'FontWeight','bold'); CH = colorbar; CH.Label.String = 'rad'; axis square; set(gca,'FontSize',fontSize,'FontWeight','normal');
set(figureHandleC,'units','normalized','position',[0,0,1,1]);