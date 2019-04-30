function [closestValues,conditionalIndices] = f_findClosestValues(valuesArray,centralValue,proximity)
differences = abs(centralValue-valuesArray);
sortedDifferences = sort(differences);
conditionalIndices = ismember(differences,sortedDifferences(1:proximity));
closestValues = valuesArray(conditionalIndices);
end