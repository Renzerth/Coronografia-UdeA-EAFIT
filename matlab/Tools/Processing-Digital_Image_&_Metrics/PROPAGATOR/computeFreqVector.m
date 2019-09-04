function [freqVectX, freqVectY,spatialXPeriod,spatialYPeriod,analysisScaling,normNMFactor,synthesisScaling] = computeFreqVector(planeSize,spatialSamples, samplingFactor)
%% Compute frequency properties
% spanFactor = 2; % Domain representation size -> Twice for padding (better sampling), once for equal size (less resolution)
freqXsamples = samplingFactor*spatialSamples(1); % Frequency domain x sampling bins
freqYsamples = samplingFactor*spatialSamples(2); % Frequency domain y sampling bins

spatialXPeriod = planeSize(1)/spatialSamples(1);  % x spatial period sampling
spatialYPeriod = planeSize(2)/spatialSamples(2);  % y spatial period sampling

samplingFreqX = 1/spatialXPeriod; % sampling frequency of x space domain - Signal's Frequency Sampling rate
samplingFreqY = 1/spatialYPeriod; % sampling frequency of y space domain - Signal's Frequency Sampling rate

DFX = samplingFreqX/(freqXsamples); % Frequency domain x components resolution
DFY = samplingFreqY/(freqYsamples); % Frequency domain y components resolution
%% Compute frequency vector
freqVectX = (-samplingFreqX/2:DFX:samplingFreqX/2-DFX) + mod(freqXsamples,2)*DFX/2; % Nyquist Range for x frequencies
freqVectY = (-samplingFreqY/2:DFY:samplingFreqY/2-DFY) + mod(freqYsamples,2)*DFY/2; % Nyquist Range for y frequencies
%% DFT Scaling parameters
analysisScaling = spatialXPeriod*spatialYPeriod;
normNMFactor = numel(spatialSamples);
synthesisScaling = 1/(analysisScaling*normNMFactor);
end