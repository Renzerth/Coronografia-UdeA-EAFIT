function [dataCoordX, dataCoordY, freqMeshX, freqMeshY, freqRadii, freqCirc] = prepareSpectralProp_mod(spaceSize, freqVectX, freqVectY, lambda)
%% Padding Window
freqXsamples = 2*spaceSize(1);
freqYsamples = 2*spaceSize(2);

addOddX = mod(spaceSize(1),2)/2;
addOddY = mod(spaceSize(2),2)/2;

dataCoordX = (freqXsamples/4+1-addOddX):(freqXsamples*3/4-addOddX);
dataCoordY = (freqYsamples/4+1-addOddY):(freqYsamples*3/4-addOddY);

%% Frequency Mesh Coordinates
[freqMeshX, freqMeshY] = meshgrid(freqVectX, freqVectY);
freqRadialFreq = freqMeshX.^2 + freqMeshY.^2;
freqRadii = sqrt(freqRadialFreq);
freqCirc = freqRadii < 1/lambda; % Evanecent waves restriction
end