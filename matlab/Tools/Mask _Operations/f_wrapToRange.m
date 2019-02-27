function [dataArray] = f_wrapToRange(dataArray, minValue, maxValue)
% Wraps a  data array in a specific range given by min and max Value
if minValue >= maxValue
    error('Value ranges are not valid. Either are equal or minValue is larger than maxValue.')
end

wrapRange = maxValue - minValue;
inRange = (dataArray < minValue) | (maxValue < dataArray);
dataArray(inRange) = mod(dataArray(inRange) - minValue, wrapRange) + minValue;
end