function [scalingLimits] = computeScaleRange(referentialData, scaleFactor, scaleShift)
finiteDataValues = referentialData(~isinf(referentialData));
refMax = scaleFactor(1)*max(finiteDataValues) + scaleShift(1);
refMin = scaleFactor(2)*min(finiteDataValues) + scaleShift(2);
scalingLimits = [refMin, refMax];
end