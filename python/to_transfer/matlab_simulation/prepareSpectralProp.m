function [paddedSpectrum, dataCoordX, dataCoordY, freqMeshX, freqMeshY, freqCirc] = prepareSpectralProp(spaceSize,samplingFactor, inputPlane, freqVectX, freqVectY, analysisScaling, lambda)
%% Padding Window
freqXsamples = samplingFactor*spaceSize(1);
freqYsamples = samplingFactor*spaceSize(2);
paddingWindow = zeros(freqXsamples,freqYsamples);
if samplingFactor == 2
    addOddX = mod(spaceSize(1),2)/2;
    addOddY = mod(spaceSize(2),2)/2;
    dataCoordX = (freqXsamples/4+1-addOddX):(freqXsamples*3/4-addOddX);
    dataCoordY = (freqYsamples/4+1-addOddY):(freqYsamples*3/4-addOddY);
else
    dataCoordX = linspace(1,freqXsamples,freqXsamples);
    dataCoordY = linspace(1,freqYsamples,freqYsamples);
end
paddingWindow(dataCoordX,dataCoordY) = inputPlane;
paddedSpectrum= fft2s(paddingWindow,analysisScaling);

%% Frequency Mesh Coordinates
[freqMeshX, freqMeshY] = meshgrid(freqVectX, freqVectY);
freqRadialFreq = freqMeshX.^2 + freqMeshY.^2;
freqRadii = sqrt(freqRadialFreq);
freqCirc = freqRadii < 1/lambda;
end
