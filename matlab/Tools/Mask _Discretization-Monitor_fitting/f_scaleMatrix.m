% Rescales each matrix data on the range [scaleMin,scaleMax]
function [rescaledData] = f_scaleMatrix(matrixData,scaleMin,scaleMax)
minValue = min(matrixData(:));
maxValue = max(matrixData(:));

rescaledData = scaleMin + (matrixData-minValue)./...
              (maxValue-minValue).*(scaleMax - scaleMin);
end