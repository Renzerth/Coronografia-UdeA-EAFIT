function [cropRange] = f_computePSFCropRange(rangeFactor,diskSizePx,centerPoint)
% diskSizePx: diameter
%% Compute visualization range for PSF intensity
cropFactor = rangeFactor*diskSizePx;
rowRange = (centerPoint(2) - cropFactor):(centerPoint(2) + cropFactor);
colRange = (centerPoint(1) - cropFactor):(centerPoint(1) + cropFactor);
cropRange = ceil([rowRange; colRange]);
end
croppedData = PSFdata(cropRange(2,:),cropRange(1,:));
rangeFactor,aproxRadius,aproxCenter