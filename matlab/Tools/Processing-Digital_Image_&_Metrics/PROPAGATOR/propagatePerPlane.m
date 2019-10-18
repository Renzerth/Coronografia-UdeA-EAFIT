function [propagationIntensities,transverseIntensityCut] = propagatePerPlane(targetPlane, k, angularSpectrum, propagationKernel, spaceSamples,halfSize,transverseCutRange,cutDataLength,propagationDistances,totalPlanes,normNMFactor, analysisScaling, synthesisScaling, dataCoordX, dataCoordY, freqCirc)
%% Initializing
propagationIntensities = zeros(spaceSamples(1),spaceSamples(2),totalPlanes);
transverseIntensityCut = zeros(cutDataLength,totalPlanes);
%% Propagate per distance
for distanceIndex = 1: totalPlanes
    localKernel = propagationKernel(freqCirc, k, angularSpectrum, propagationDistances(distanceIndex));
    propagationIntensities(:,:,distanceIndex) = convoluteSignal(localKernel, targetPlane, normNMFactor, analysisScaling, synthesisScaling, dataCoordX, dataCoordY, halfSize);
    transverseIntensityCut(:,distanceIndex) = abs(propagationIntensities(halfSize(2),transverseCutRange,distanceIndex)).^2;
end
end