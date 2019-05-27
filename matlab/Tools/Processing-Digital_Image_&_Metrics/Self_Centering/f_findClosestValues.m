function [closestValues,conditionalIndices] = f_findClosestValues(valuesArray,centralValue,proximity)
% Finds the closest values inside a vector within the first N approximated
% values
differences = abs(centralValue-valuesArray);
sortedDifferences = sort(differences);
conditionalIndices = ismember(differences,sortedDifferences(1:proximity));
closestValues = valuesArray(conditionalIndices);
end