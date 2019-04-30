function [multiple] = f_getClosestMultiple(number,multiple)

rangeVals = 1:ceil(number/2);
multiples = [rangeVals(rem(number,rangeVals)==0), number];
[multiple,~] = min(f_findClosestValues(multiples,multiple,1));
% multiple = multiple - rem(number,multiple);
% if rem(number,multiple) ~= 0
%     multiple = multiple - round(rem(number,multiple/2));
% end
end