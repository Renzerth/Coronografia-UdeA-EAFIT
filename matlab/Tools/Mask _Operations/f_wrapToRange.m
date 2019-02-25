function [dataArray] = f_wrapToRange(dataArray, minValue, maxValue)
if minValue >= maxValue
    error('Value ranges are not valid. Either are equal or minValue is larger than maxValue.')
end

wrapRange = maxValue - minValue;
inRange = (dataArray < minValue) | (maxValue < dataArray);
dataArray(inRange) = mod(dataArray(inRange) - minValue, wrapRange) + minValue;
end