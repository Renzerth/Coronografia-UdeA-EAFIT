function [dataArray] = f_wrapToRange(dataArray, minValue, maxValue)
% Wraps a  data array in a specific range given by min and max Value
% It doesn't work for wrapped structures that are inputs, meaning that the
% dataArray should have values bigger than the min,max-Value ranges
% Inputs:
%  dataArray: n-dimensional numerical array that should have a range
%  outside the [minValue, maxValue] range, otherwise the data wouldn't wrap
%  [minValue, maxValue]: range for the wrapping
%
% Output:
%  wrapped data
%
if minValue >= maxValue
    error('Value ranges are not valid. Either are equal or minValue is larger than maxValue.')
end

wrapRange = maxValue - minValue;
inRange = (dataArray < minValue) | (maxValue < dataArray);
dataArray(inRange) = mod(dataArray(inRange) - minValue, wrapRange) + minValue;
end