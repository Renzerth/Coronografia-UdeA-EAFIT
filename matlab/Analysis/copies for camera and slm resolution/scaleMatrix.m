function [rescaledData] = scaleMatrix(matrixData,scaleMin,scaleMax)
minValue = min(matrixData(:));
maxValue = max(matrixData(:));

rescaledData = scaleMin + (matrixData-minValue)./(maxValue-minValue).*(scaleMax - scaleMin);
end