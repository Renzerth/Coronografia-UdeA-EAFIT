function [xCoordinates, yCoordinates, X, Y, Theta, Rho] = makeSpaceCoordinates(xSize,ySize,xSamples,ySamples,midPointX,midPointY)

samplingXPeriod = xSize/xSamples; % space sampling period
samplingYPeriod = ySize/ySamples; % space sampling period

DX = 1; % Digital sampling period
DY = 1; % Digital sampling period

halfXSamples = floor((xSamples + 1)/2) - mod(xSamples,2);
halfYSamples = floor((ySamples + 1)/2) - mod(ySamples,2);

centerShiftX = (midPointX - halfXSamples);
centerShiftY = (midPointY - halfYSamples);

xCoordinates = samplingXPeriod*((-halfXSamples:DX:halfXSamples - DX + mod(xSamples,2)*DX) - (centerShiftX));
yCoordinates = samplingYPeriod*((-halfYSamples:DY:halfYSamples - DY + mod(ySamples,2)*DY) - (centerShiftY));

[X,Y] = meshgrid(xCoordinates, yCoordinates);
[Theta, Rho] = cart2pol(X,Y);

end