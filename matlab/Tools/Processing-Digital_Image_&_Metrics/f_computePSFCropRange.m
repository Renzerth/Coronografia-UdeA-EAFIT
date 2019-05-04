function [cropRange] = f_computePSFCropRange(rangeFactor,diskSizePx,centerPoint)
%% Compute visualization range for PSF intensity
cropFactor = rangeFactor*diskSizePx;
rowRange = (centerPoint(2) - cropFactor):(centerPoint(2) + cropFactor);
colRange = (centerPoint(1) - cropFactor):(centerPoint(1) + cropFactor);
cropRange = [rowRange; colRange];
end